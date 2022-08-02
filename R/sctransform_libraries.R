# @useDynLib sctransform
#NULL


#' Variance stabilizing transformation for UMI count data
#'
#' Apply variance stabilizing transformation to UMI count data using a regularized Negative Binomial regression model.
#' This will remove unwanted effects from UMI data and return Pearson residuals.
#' Uses future_lapply; you can set the number of cores it will use to n with plan(strategy = "multicore", workers = n).
#' If n_genes is set, only a (somewhat-random) subset of genes is used for estimating the
#' initial model parameters. For details see \href{http://dx.doi.org/10.1186/s13059-019-1874-1}{doi: 10.1186/s13059-019-1874-1}.
#'
#' @param umi A matrix of UMI counts with genes as rows and cells as columns
#' @param cell_attr A data frame containing the dependent variables; if omitted a data frame with umi and gene will be generated
#' @param latent_var The independent variables to regress out as a character vector; must match column names in cell_attr; default is c("log_umi")
#' @param batch_var The dependent variables indicating which batch a cell belongs to; no batch interaction terms used if omiited
#' @param latent_var_nonreg The non-regularized dependent variables to regress out as a character vector; must match column names in cell_attr; default is NULL
#' @param n_genes Number of genes to use when estimating parameters (default uses 2000 genes, set to NULL to use all genes)
#' @param n_cells Number of cells to use when estimating parameters (default uses all cells)
#' @param method Method to use for initial parameter estimation; one of 'poisson', 'qpoisson', 'nb_fast', 'nb', 'nb_theta_given', 'glmGamPoi', 'offset', 'offset_shared_theta_estimate'; default is 'poisson'
#' @param do_regularize Boolean that, if set to FALSE, will bypass parameter regularization and use all genes in first step (ignoring n_genes); default is FALSE
#' @param theta_regularization Method to use to regularize theta; use 'log_theta' for the behavior prior to version 0.3; default is 'od_factor'
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param bin_size Number of genes to process simultaneously; this will determine how often the progress bars are updated and how much memory is being used; default is 500
#' @param min_cells Only use genes that have been detected in at least this many cells; default is 5
#' @param residual_type What type of residuals to return; can be 'pearson', 'deviance', or 'none'; default is 'pearson'
#' @param return_cell_attr Make cell attributes part of the output; default is FALSE
#' @param return_gene_attr Calculate gene attributes and make part of output; default is TRUE
#' @param return_corrected_umi If set to TRUE output will contain corrected UMI matrix; see \code{correct} function
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is -Inf
#' @param bw_adjust Kernel bandwidth adjustment factor used during regurlarization; factor will be applied to output of bw.SJ; default is 3
#' @param gmean_eps Small value added when calculating geometric mean of a gene to avoid log(0); default is 1
#' @param theta_estimation_fun Character string indicating which method to use to estimate theta (when method = poisson); default is 'theta.ml', but 'theta.mm' seems to be a good and fast alternative
#' @param theta_given If method is set to nb_theta_given, this should be a named numeric vector of fixed theta values for the genes; if method is offset, this should be a single value; default is NULL
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A list with components
#' \item{y}{Matrix of transformed data, i.e. Pearson residuals, or deviance residuals; empty if \code{residual_type = 'none'}}
#' \item{umi_corrected}{Matrix of corrected UMI counts (optional)}
#' \item{model_str}{Character representation of the model formula}
#' \item{model_pars}{Matrix of estimated model parameters per gene (theta and regression coefficients)}
#' \item{model_pars_outliers}{Vector indicating whether a gene was considered to be an outlier}
#' \item{model_pars_fit}{Matrix of fitted / regularized model parameters}
#' \item{model_str_nonreg}{Character representation of model for non-regularized variables}
#' \item{model_pars_nonreg}{Model parameters for non-regularized variables}
#' \item{genes_log_gmean_step1}{log-geometric mean of genes used in initial step of parameter estimation}
#' \item{cells_step1}{Cells used in initial step of parameter estimation}
#' \item{arguments}{List of function call arguments}
#' \item{cell_attr}{Data frame of cell meta data (optional)}
#' \item{gene_attr}{Data frame with gene attributes such as mean, detection rate, etc. (optional)}
#' \item{times}{Time stamps at various points in the function}
#'
#' @section Details:
#' In the first step of the algorithm, per-gene glm model parameters are learned. This step can be done
#' on a subset of genes and/or cells to speed things up.
#' If \code{method} is set to 'poisson', a poisson regression is done and
#' the negative binomial theta parameter is estimated using the response residuals in
#' \code{theta_estimation_fun}.
#' If \code{method} is set to 'qpoisson', coefficients and overdispersion (phi) are estimated by quasi 
#' poisson regression and theta is estimated based on phi and the mean fitted value - this is currently 
#' the fastest method with results very similar to 'glmGamPoi'
#' If \code{method} is set to 'nb_fast', coefficients and theta are estimated as in the
#' 'poisson' method, but coefficients are then re-estimated using a proper negative binomial
#' model in a second call to glm with \code{family = MASS::negative.binomial(theta = theta)}.
#' If \code{method} is set to 'nb', coefficients and theta are estimated by a single call to
#' \code{MASS::glm.nb}.
#' If \code{method} is set to 'glmGamPoi', coefficients and theta are estimated by a single call to
#' \code{glmGamPoi::glm_gp}.
#' 
#' A special case is \code{method = 'offset'}. Here no regression parameters are learned, but
#' instead an offset model is assumed. The latent variable is set to log_umi and a fixed 
#' slope of log(10) is used (offset). The intercept is given by log(gene_mean) - log(avg_cell_umi). 
#' See Lause et al. (\href{https://doi.org/10.1101/2020.12.01.405886}{bioRxiv 2020.12.01.405886}) for details.
#' Theta is set
#' to 100 by default, but can be changed using the \code{theta_given} parameter (single numeric value).
#' If the offset method is used, the following parameters are overwritten:
#' \code{cell_attr <- NULL, latent_var <- c('log_umi'), batch_var <- NULL, latent_var_nonreg <- NULL,
#' n_genes <- NULL, n_cells <- NULL, do_regularize <- FALSE}. Further, \code{method = 'offset_shared_theta_estimate'}
#' exists where the 250 most highly expressed genes with detection rate of at least 0.5 are used
#' to estimate a theta that is then shared across all genes. Thetas are estimated per individual gene
#' using 5000 randomly selected cells. The final theta used for all genes is then the average.
#' 
#'
#' @import Matrix
#' @importFrom future.apply future_lapply
#' @importFrom MASS theta.ml theta.mm glm.nb negative.binomial
#' @importFrom stats glm glm.fit df.residual ksmooth model.matrix as.formula approx density poisson var bw.SJ
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc)
#' }
#'
vst <- function(umi,
                cell_attr = NULL,
                latent_var = c('log_umi'),
                batch_var = NULL,
                latent_var_nonreg = NULL,
                n_genes = 2000,
                n_cells = NULL,
                method = 'poisson',
                do_regularize = TRUE,
                theta_regularization = 'od_factor',
                res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                bin_size = 500,
                min_cells = 5,
                residual_type = 'pearson',
                return_cell_attr = FALSE,
                return_gene_attr = TRUE,
                return_corrected_umi = FALSE,
                min_variance = -Inf,
                bw_adjust = 3,
                gmean_eps = 1,
                theta_estimation_fun = 'theta.ml',
                theta_given = NULL,
                verbosity = 2,
                verbose = NULL,
                show_progress = NULL) {
  arguments <- as.list(environment())
  arguments <- arguments[!names(arguments) %in% c("umi", "cell_attr")]
  
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }
  
  # Check for suggested package
  if (method == "glmGamPoi") {
    glmGamPoi_check <- requireNamespace("glmGamPoi", quietly = TRUE)
    if (!glmGamPoi_check){
      stop('Please install the glmGamPoi package. See https://github.com/const-ae/glmGamPoi for details.')
    }
  }
  
  # Special case offset model - override most parameters
  if (startsWith(x = method, prefix = 'offset')) {
    cell_attr <- NULL
    latent_var <- c('log_umi')
    batch_var <- NULL
    latent_var_nonreg <- NULL
    n_genes <- NULL
    n_cells <- NULL
    do_regularize <- FALSE
    if (is.null(theta_given)) {
      theta_given <- 100
    } else {
      theta_given <- theta_given[1]
    }
  }
  
  
  times <- list(start_time = Sys.time())
  
  cell_attr <- make_cell_attr(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity)
  if (!is.null(batch_var)) {
    cell_attr[, batch_var] <- as.factor(cell_attr[, batch_var])
    batch_levels <- levels(cell_attr[, batch_var])
  }
  
  # we will generate output for all genes detected in at least min_cells cells
  # but for the first step of parameter estimation we might use only a subset of genes
  genes_cell_count <- rowSums(umi >= 0.01)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))
  
  if (!do_regularize && !is.null(n_genes)) {
    if (verbosity > 0) {
      message('do_regularize is set to FALSE, will use all genes')
    }
    n_genes <- NULL
  }
  
  if (!is.null(n_cells) && n_cells < ncol(umi)) {
    # downsample cells to speed up the first step
    cells_step1 <- sample(x = colnames(umi), size = n_cells)
    if (!is.null(batch_var)) {
      dropped_batch_levels <- setdiff(batch_levels, levels(droplevels(cell_attr[cells_step1, batch_var])))
      if (length(dropped_batch_levels) > 0) {
        stop('Dropped batch levels ', dropped_batch_levels, ', set n_cells higher')
      }
    }
    genes_cell_count_step1 <- rowSums(umi[, cells_step1] > 0)
    genes_step1 <- rownames(umi)[genes_cell_count_step1 >= min_cells]
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  } else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean
  }
  
  data_step1 <- cell_attr[cells_step1, , drop = FALSE]
  
  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    # density-sample genes to speed up the first step
    log_gmean_dens <- density(x = genes_log_gmean_step1, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_gmean_dens$x, y = log_gmean_dens$y, xout = genes_log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, prob = sampling_prob)
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  }
  
  if (!is.null(batch_var)) {
    model_str <- paste0('y ~ (', paste(latent_var, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
  } else {
    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))
  }
  
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Variance stabilizing transformation of count matrix of size ', nrow(umi), ' by ', ncol(umi))
    message('Model formula is ', model_str)
  }
  
  times$get_model_pars = Sys.time()
  model_pars <- get_model_pars(genes_step1, bin_size, umi, model_str, cells_step1,
                               method, data_step1, theta_given, theta_estimation_fun,
                               verbosity)
  # make sure theta is not too small
  min_theta <- 1e-7
  if (any(model_pars[, 'theta'] < min_theta)) {
    if (verbosity > 0) {
      msg <- sprintf('There are %d estimated thetas smaller than %g - will be set to %g', sum(model_pars[, 'theta'] < min_theta), min_theta, min_theta)
      message(msg)
    }
    model_pars[, 'theta'] <- pmax(model_pars[, 'theta'], min_theta)
  }
  
  
  times$reg_model_pars = Sys.time()
  if (do_regularize) {
    model_pars_fit <- reg_model_pars(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                                     batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                                     theta_regularization, verbosity)
    model_pars_outliers <- attr(model_pars_fit, 'outliers')
  } else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }
  
  # use all fitted values in NB model
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)
  
  if (!is.null(latent_var_nonreg)) {
    if (verbosity > 0) {
      message('Estimating parameters for following non-regularized variables: ', latent_var_nonreg)
    }
    if (!is.null(batch_var)) {
      model_str_nonreg <- paste0('y ~ (', paste(latent_var_nonreg, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
    } else {
      model_str_nonreg <- paste0('y ~ ', paste(latent_var_nonreg, collapse = ' + '))
    }
    
    times$get_model_pars_nonreg = Sys.time()
    model_pars_nonreg <- get_model_pars_nonreg(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, verbosity)
    
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', model_str_nonreg)), cell_attr)
    model_pars_final <- cbind(model_pars_fit, model_pars_nonreg)
    regressor_data_final <- cbind(regressor_data, regressor_data_nonreg)
    #model_pars_final[, '(Intercept)'] <- model_pars_final[, '(Intercept)'] + model_pars_nonreg[, '(Intercept)']
    #model_pars_final <- cbind(model_pars_final, model_pars_nonreg[, -1, drop=FALSE])
    # model_str <- paste0(model_str, gsub('^y ~ 1', '', model_str2))
  } else {
    model_str_nonreg <- ''
    model_pars_nonreg <- c()
    model_pars_final <- model_pars_fit
    regressor_data_final <- regressor_data
  }
  
  times$get_residuals = Sys.time()
  if (!residual_type == 'none') {
    if (verbosity > 0) {
      message('Second step: Get residuals using fitted parameters for ', length(x = genes), ' genes')
    }
    bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
    max_bin <- max(bin_ind)
    if (verbosity > 1) {
      pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    }
    res <- matrix(NA_real_, length(genes), nrow(regressor_data_final), dimnames = list(genes, rownames(regressor_data_final)))
    for (i in 1:max_bin) {
      genes_bin <- genes[bin_ind == i]
      mu <- exp(tcrossprod(model_pars_final[genes_bin, -1, drop=FALSE], regressor_data_final))
      y <- as.matrix(umi[genes_bin, , drop=FALSE])
      res[genes_bin, ] <- switch(residual_type,
                                 'pearson' = pearson_residual(y, mu, model_pars_final[genes_bin, 'theta'], min_var = min_variance),
                                 'deviance' = deviance_residual(y, mu, model_pars_final[genes_bin, 'theta']),
                                 stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment')
      )
      if (verbosity > 1) {
        setTxtProgressBar(pb, i)
      }
    }
    if (verbosity > 1) {
      close(pb)
    }
  } else {
    if (verbosity > 0) {
      message('Skip calculation of full residual matrix')
    }
    res <- matrix(data = NA, nrow = 0, ncol = 0)
  }
  
  rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_outliers = model_pars_outliers,
             model_pars_fit = model_pars_fit,
             model_str_nonreg = model_str_nonreg,
             model_pars_nonreg = model_pars_nonreg,
             arguments = arguments,
             genes_log_gmean_step1 = genes_log_gmean_step1,
             cells_step1 = cells_step1,
             cell_attr = cell_attr)
  rm(res)
  gc(verbose = FALSE)
  
  times$correct_umi = Sys.time()
  if (return_corrected_umi) {
    if (residual_type != 'pearson') {
      message("Will not return corrected UMI because residual type is not set to 'pearson'")
    } else {
      rv$umi_corrected <- sctransform::correct(rv, do_round = TRUE, do_pos = TRUE,
                                               verbosity = verbosity)
      rv$umi_corrected <- as(object = rv$umi_corrected, Class = 'dgCMatrix')
    }
  }
  
  rv$y[rv$y < res_clip_range[1]] <- res_clip_range[1]
  rv$y[rv$y > res_clip_range[2]] <- res_clip_range[2]
  
  if (!return_cell_attr) {
    rv[['cell_attr']] <- NULL
  }
  
  times$get_gene_attr = Sys.time()
  if (return_gene_attr) {
    if (verbosity > 0) {
      message('Calculating gene attributes')
    }
    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      gmean = 10 ^ genes_log_gmean,
      variance = row_var(umi))
    if (ncol(rv$y) > 0) {
      gene_attr$residual_mean = rowMeans(rv$y)
      gene_attr$residual_variance = row_var(rv$y)
    }
    # Special case offset model - also calculate arithmetic mean
    if (startsWith(x = method, prefix = 'offset')) {
      gene_attr$amean <- rowMeans(umi)
    }
    
    rv[['gene_attr']] <- gene_attr
  }
  
  if (verbosity > 0) {
    message('Wall clock passed: ', capture.output(print(Sys.time() - times$start_time)))
  }
  times$done = Sys.time()
  rv$times <- times
  return(rv)
}


get_model_pars <- function(genes_step1, bin_size, umi, model_str, cells_step1,
                           method, data_step1, theta_given, theta_estimation_fun,
                           verbosity) {
  # Special case offset model with one theta for all genes
  if (startsWith(x = method, prefix = 'offset')) {
    gene_mean <- rowMeans(umi)
    mean_cell_sum <- mean(colSums(umi))
    model_pars <- cbind(rep(theta_given, nrow(umi)),
                        log(gene_mean) - log(mean_cell_sum),
                        rep(log(10), nrow(umi)))
    dimnames(model_pars) <- list(rownames(umi), c('theta', '(Intercept)', 'log_umi'))
    if (method == 'offset_shared_theta_estimate') {
      # use all genes with detection rate > 0.5 to estimate theta
      # if there are more, use the 250 most highly expressed ones
      # use at most 5000 cells (random sample)
      use_genes <- rowMeans(umi > 0) > 0.5
      if (sum(use_genes) > 250) {
        o <- order(-row_gmean(umi[use_genes, ]))
        use_genes <- which(use_genes)[o[1:250]]
      }
      use_cells <- sample(x = ncol(umi), size = min(ncol(umi), 5000), replace = FALSE)
      if (verbosity > 0) {
        message(sprintf('Estimate shared theta for offset model using %d genes, %d cells', 
                        length(x = use_genes), length(x = use_cells)))
      }
      y <- as.matrix(umi[use_genes, use_cells])
      regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data_step1[use_cells, ])
      mu <- exp(tcrossprod(model_pars[use_genes, -1, drop=FALSE], regressor_data))
      theta <- sapply(1:nrow(y), function(i) {
        as.numeric(MASS::theta.ml(y = y[i, ], mu = mu[i, ], limit = 100))
      })
      model_pars[, 'theta'] <- mean(theta)
    }
    return(model_pars)
  }
  
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Get Negative Binomial regression parameters per gene')
    message('Using ', length(x = genes_step1), ' genes, ', length(x = cells_step1), ' cells')
  }
  
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars <- list()
  for (i in 1:max_bin) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, cells_step1, drop=FALSE])
    if (!is.null(theta_given)) {
      theta_given_bin <- theta_given[genes_bin_regress]
    }
    
    # umi_bin is a matrix of counts - we want a model per row
    # if there are multiple workers, split up the matrix in chunks of n rows
    # where n is the number of workers
    n_workers <- 1
    if (future::supportsMulticore()) {
      n_workers <- future::nbrOfWorkers()
    }
    genes_per_worker <- nrow(umi_bin) / n_workers + .Machine$double.eps
    index_vec <- 1:nrow(umi_bin)
    index_lst <- split(index_vec, ceiling(index_vec/genes_per_worker))
    
    # the index list will have at most n_workers entries, each one defining which genes to work on
    par_lst <- future_lapply(
      X = index_lst,
      FUN = function(indices) {
        umi_bin_worker <- umi_bin[indices, , drop = FALSE]
        if (method == 'poisson') {
          return(fit_poisson(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_estimation_fun = theta_estimation_fun))
        }
        if (method == 'qpoisson') {
          return(fit_qpoisson(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
        if (method == 'nb_theta_given') {
          theta_given_bin_worker <- theta_given_bin[indices]
          return(fit_nb_theta_given(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_given = theta_given_bin_worker))
        }
        if (method == 'nb_fast') {
          return(fit_nb_fast(umi = umi_bin_worker, model_str = model_str, data = data_step1, theta_estimation_fun = theta_estimation_fun))
        }
        if (method == 'nb') {
          return(fit_nb(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
        if (method == "glmGamPoi") {
          return(fit_glmGamPoi(umi = umi_bin_worker, model_str = model_str, data = data_step1))
        }
      }
    )
    model_pars[[i]] <- do.call(rbind, par_lst)
    
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  model_pars <- do.call(rbind, model_pars)
  if (verbosity > 1) {
    close(pb)
  }
  rownames(model_pars) <- genes_step1
  colnames(model_pars)[1] <- 'theta'
  return(model_pars)
}

get_model_pars_nonreg <- function(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, verbosity) {
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  model_pars_nonreg <- list()
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- tcrossprod(model_pars_fit[genes_bin, -1, drop=FALSE], regressor_data)
    umi_bin <- as.matrix(umi[genes_bin, ])
    model_pars_nonreg[[i]] <- do.call(rbind,
                                      future_lapply(genes_bin, function(gene) {
                                        fam <- negative.binomial(theta = model_pars_fit[gene, 'theta'], link = 'log')
                                        y <- umi_bin[gene, ]
                                        offs <- mu[gene, ]
                                        fit <- glm(as.formula(model_str_nonreg), data = cell_attr, family = fam, offset=offs)
                                        return(fit$coefficients)
                                      }))
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  model_pars_nonreg <- do.call(rbind, model_pars_nonreg)
  rownames(model_pars_nonreg) <- genes
  return(model_pars_nonreg)
}

reg_model_pars <- function(model_pars, genes_log_gmean_step1, genes_log_gmean, cell_attr,
                           batch_var, cells_step1, genes_step1, umi, bw_adjust, gmean_eps,
                           theta_regularization, verbosity) {
  genes <- names(genes_log_gmean)
  
  # we don't regularize theta directly
  # prior to v0.3 we regularized log10(theta)
  # now we transform to overdispersion factor
  # variance of NB is mu * (1 + mu / theta)
  # (1 + mu / theta) is what we call overdispersion factor here
  dispersion_par <- switch(theta_regularization,
                           'log_theta' = log10(model_pars[, 'theta']),
                           'od_factor' = log10(1 + 10^genes_log_gmean_step1 / model_pars[, 'theta']),
                           stop('theta_regularization ', theta_regularization, ' unknown - only log_theta and od_factor supported at the moment')
  )
  
  model_pars <- model_pars[, colnames(model_pars) != 'theta']
  model_pars <- cbind(dispersion_par, model_pars)
  
  # look for outliers in the parameters
  # outliers are those that do not fit the overall relationship with the mean at all
  outliers <- apply(model_pars, 2, function(y) is_outlier(y, genes_log_gmean_step1))
  outliers <- apply(outliers, 1, any)
  if (sum(outliers) > 0) {
    if (verbosity > 0) {
      message('Found ', sum(outliers), ' outliers - those will be ignored in fitting/regularization step\n')
    }
    model_pars <- model_pars[!outliers, ]
    genes_step1 <- rownames(model_pars)
    genes_log_gmean_step1 <- genes_log_gmean_step1[!outliers]
  }
  
  # select bandwidth to be used for smoothing
  bw <- bw.SJ(genes_log_gmean_step1) * bw_adjust
  
  # for parameter predictions
  x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1))
  x_points <- pmin(x_points, max(genes_log_gmean_step1))
  
  # take results from step 1 and fit/predict parameters to all genes
  o <- order(x_points)
  model_pars_fit <- matrix(NA_real_, length(genes), ncol(model_pars),
                           dimnames = list(genes, colnames(model_pars)))
  
  # fit / regularize dispersion parameter
  model_pars_fit[o, 'dispersion_par'] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, 'dispersion_par'],
                                                 x.points = x_points, bandwidth = bw, kernel='normal')$y
  
  if (is.null(batch_var)){
    # global fit / regularization for all coefficients
    for (i in 2:ncol(model_pars)) {
      model_pars_fit[o, i] <- ksmooth(x = genes_log_gmean_step1, y = model_pars[, i],
                                      x.points = x_points, bandwidth = bw, kernel='normal')$y
    }
  } else {
    # fit / regularize per batch
    batches <- unique(cell_attr[, batch_var])
    for (b in batches) {
      sel <- cell_attr[, batch_var] == b & rownames(cell_attr) %in% cells_step1
      #batch_genes_log_gmean_step1 <- log10(rowMeans(umi[genes_step1, sel]))
      batch_genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, sel], eps = gmean_eps))
      if (any(is.infinite(batch_genes_log_gmean_step1))) {
        if (verbosity > 0) {
          message('Some genes not detected in batch ', b, ' -- assuming a low mean.')
        }
        batch_genes_log_gmean_step1[is.infinite(batch_genes_log_gmean_step1) & batch_genes_log_gmean_step1 < 0] <- min(batch_genes_log_gmean_step1[!is.infinite(batch_genes_log_gmean_step1)])
      }
      sel <- cell_attr[, batch_var] == b
      #batch_genes_log_gmean <- log10(rowMeans(umi[, sel]))
      batch_genes_log_gmean <- log10(row_gmean(umi[, sel], eps = gmean_eps))
      # in case some genes have not been observed in this batch
      batch_genes_log_gmean <- pmax(batch_genes_log_gmean, min(batch_genes_log_gmean_step1))
      batch_o <- order(batch_genes_log_gmean)
      for (i in which(grepl(paste0(batch_var, b), colnames(model_pars)))) {
        model_pars_fit[batch_o, i] <- ksmooth(x = batch_genes_log_gmean_step1, y = model_pars[, i],
                                              x.points = batch_genes_log_gmean, bandwidth = bw, kernel='normal')$y
      }
    }
  }
  
  # back-transform dispersion parameter to theta
  theta <- switch(theta_regularization,
                  'log_theta' = 10^model_pars_fit[, 'dispersion_par'],
                  'od_factor' = 10^genes_log_gmean / (10^model_pars_fit[, 'dispersion_par'] - 1)
  )
  model_pars_fit <- model_pars_fit[, colnames(model_pars_fit) != 'dispersion_par']
  model_pars_fit <- cbind(theta, model_pars_fit)
  
  attr(model_pars_fit, 'outliers') <- outliers
  return(model_pars_fit)
}


# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

row_mean_dgcmatrix <- function(matrix) {
  .Call('_sctransform_row_mean_dgcmatrix', PACKAGE = 'sctransform', matrix)
}

row_mean_grouped_dgcmatrix <- function(matrix, group, shuffle) {
  .Call('_sctransform_row_mean_grouped_dgcmatrix', PACKAGE = 'sctransform', matrix, group, shuffle)
}

row_gmean_dgcmatrix <- function(matrix, eps) {
  .Call('_sctransform_row_gmean_dgcmatrix', PACKAGE = 'sctransform', matrix, eps)
}

row_gmean_grouped_dgcmatrix <- function(matrix, group, eps, shuffle) {
  .Call('_sctransform_row_gmean_grouped_dgcmatrix', PACKAGE = 'sctransform', matrix, group, eps, shuffle)
}

row_nonzero_count_dgcmatrix <- function(matrix) {
  .Call('_sctransform_row_nonzero_count_dgcmatrix', PACKAGE = 'sctransform', matrix)
}

row_nonzero_count_grouped_dgcmatrix <- function(matrix, group) {
  .Call('_sctransform_row_nonzero_count_grouped_dgcmatrix', PACKAGE = 'sctransform', matrix, group)
}

row_var_dgcmatrix <- function(x, i, rows, cols) {
  .Call('_sctransform_row_var_dgcmatrix', PACKAGE = 'sctransform', x, i, rows, cols)
}

grouped_mean_diff_per_row <- function(x, group, shuffle) {
  .Call('_sctransform_grouped_mean_diff_per_row', PACKAGE = 'sctransform', x, group, shuffle)
}

mean_boot <- function(x, N, S) {
  .Call('_sctransform_mean_boot', PACKAGE = 'sctransform', x, N, S)
}

mean_boot_grouped <- function(x, group, N, S) {
  .Call('_sctransform_mean_boot_grouped', PACKAGE = 'sctransform', x, group, N, S)
}

distribution_shift <- function(x) {
  .Call('_sctransform_distribution_shift', PACKAGE = 'sctransform', x)
}

qpois_reg <- function(X, Y, tol, maxiters, minphi, returnfit) {
  .Call('_sctransform_qpois_reg', PACKAGE = 'sctransform', X, Y, tol, maxiters, minphi, returnfit)
}


#Fir NB regression models using different approaches

fit_poisson <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
                    'theta.ml' = as.numeric(x = theta.ml(y = y, mu = fit$fitted)),
                    'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
                    stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    return(c(theta, fit$coefficients))
  }))
  return(par_mat)
}

fit_qpoisson <- function(umi, model_str, data) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  par_mat <- t(apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, FALSE)
    return(c(fit$theta.guesstimate, fit$coefficients))
  }))
  return(par_mat)
}

fit_nb_theta_given <- function(umi, model_str, data, theta_given) {
  par_lst <- lapply(1:nrow(umi), function(j) {
    y <- umi[j, ]
    theta <- theta_given[j]
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, glm(as.formula(model_str), data = data, family = poisson)$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(do.call(rbind, par_lst))
}

fit_nb_fast <- function(umi, model_str, data, theta_estimation_fun) {
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), data)
  dfr <- ncol(umi) - ncol(regressor_data)
  par_mat <- apply(umi, 1, function(y) {
    fit <- qpois_reg(regressor_data, y, 1e-9, 100, 1.0001, TRUE)
    theta <- switch(theta_estimation_fun,
                    'theta.ml' = as.numeric(x = theta.ml(y = y, mu = fit$fitted)),
                    'theta.mm' = as.numeric(x = theta.mm(y = y, mu = fit$fitted, dfr = dfr)),
                    stop('theta_estimation_fun ', theta_estimation_fun, ' unknown - only theta.ml and theta.mm supported at the moment')
    )
    fit2 <- 0
    try(fit2 <- glm(as.formula(model_str), data = data, family = negative.binomial(theta=theta)), silent=TRUE)
    if (inherits(x = fit2, what = 'numeric')) {
      return(c(theta, fit$coefficients))
    } else {
      return(c(theta, fit2$coefficients))
    }
  })
  return(t(par_mat))
}

fit_nb <- function(umi, model_str, data) {
  par_mat <- apply(umi, 1, function(y) {
    fit <- 0
    try(fit <- glm.nb(as.formula(model_str), data = data), silent=TRUE)
    if (inherits(x = fit, what = 'numeric')) {
      fit <- glm(as.formula(model_str), data = data, family = poisson)
      fit$theta <- as.numeric(x = theta.ml(y = y, mu = fit$fitted))
    }
    return(c(fit$theta, fit$coefficients))
  })
  return(t(par_mat))
}

fit_glmGamPoi <- function(umi, model_str, data) {
  fit <- glmGamPoi::glm_gp(data = umi,
                           design = as.formula(gsub("y", "", model_str)),
                           col_data = data,
                           size_factors = FALSE)
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-4)
  colnames(fit$Beta)[match(x = 'Intercept', colnames(fit$Beta))] <- "(Intercept)"
  return(cbind(fit$theta, fit$Beta))
}


# Check cell attributes; add missing ones
make_cell_attr <- function(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity) {
  if (is.null(cell_attr)) {
    cell_attr <- data.frame(row.names = colnames(umi))
  }
  
  # Do not allow certain variable names
  no_good <- c('(Intercept)', 'Intercept')
  if (any(no_good %in% c(latent_var, batch_var, latent_var_nonreg))) {
    stop('Do not use the following variable names for a latent variable or batch variable: ', paste(no_good, collapse = ', '))
  }
  
  # these are the cell attributes that we know how to calculate given the count matrix
  known_attr <- c('umi', 'gene', 'log_umi', 'log_gene', 'umi_per_gene', 'log_umi_per_gene')
  # these are the missing cell attributes specified in latent_var
  missing_attr <- setdiff(c(latent_var, batch_var, latent_var_nonreg), colnames(cell_attr))
  if (length(missing_attr) > 0) {
    if (verbosity > 0) {
      message('Calculating cell attributes from input UMI matrix: ', paste(missing_attr, collapse = ', '))
    }
    unknown_attr <- setdiff(missing_attr, known_attr)
    if (length(unknown_attr) > 0) {
      stop(sprintf('Unknown cell attributes: %s. Check latent_var, batch_var and latent_var_nonreg and make sure the variables are in cell_attr', paste(unknown_attr, collapse = ', ')))
    }
    new_attr <- list()
    if (any(c('umi', 'log_umi', 'umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$umi <- colSums(umi)
      new_attr$log_umi <- log10(new_attr$umi)
    }
    if (any(c('gene', 'log_gene', 'umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$gene <- colSums(umi > 0)
      new_attr$log_gene <- log10(new_attr$gene)
    }
    if (any(c('umi_per_gene', 'log_umi_per_gene') %in% missing_attr)) {
      new_attr$umi_per_gene <- new_attr$umi / new_attr$gene
      new_attr$log_umi_per_gene <- log10(new_attr$umi_per_gene)
    }
    new_attr <- do.call(cbind, new_attr)
    cell_attr <- cbind(cell_attr, new_attr[, setdiff(colnames(new_attr), colnames(cell_attr)), drop = TRUE])
  }
  
  # make sure no NA, NaN, Inf values are in cell attributes - they would cause
  # problems later on
  rel_attr <- cell_attr[, c(latent_var, batch_var, latent_var_nonreg)]
  if (any(is.na(rel_attr)) || 
      any(is.nan(rel_attr)) || 
      any(is.infinite(rel_attr))) {
    stop('cell attributes cannot contain any NA, NaN, or infinite values')
  }
  
  return(cell_attr)
}

#' Geometric mean per row
#'
#' @param x matrix of class \code{matrix} or \code{dgCMatrix}
#' @param eps small value to add to x to avoid log(0); default is 1
#'
#' @return geometric means
row_gmean <- function(x, eps = 1) {
  if (inherits(x = x, what = 'matrix')) {
    return(exp(rowMeans(log(x + eps))) - eps)
  }
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_gmean_dgcmatrix(matrix = x, eps = eps)
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}



#' Variance per row
#'
#' @param x matrix of class \code{matrix} or \code{dgCMatrix}
#'
#' @return variances
#' 
#' @importFrom matrixStats rowVars
row_var <- function(x) {
  if (inherits(x = x, what = 'matrix')) {
    ret <- rowVars(x)
    names(ret) <- rownames(x)
    return(ret)
  }
  if (inherits(x = x, what = 'dgCMatrix')) {
    ret <- row_var_dgcmatrix(x = x@x, i = x@i, rows = nrow(x), cols = ncol(x))
    names(ret) <- rownames(x)
    return(ret)
  }
  stop('matrix x needs to be of class matrix or dgCMatrix')
}



#' Identify outliers
#'
#' @param y Dependent variable
#' @param x Independent variable
#' @param th Outlier score threshold
#'
#' @return Boolean vector
#'
#' @importFrom stats aggregate
#'
is_outlier <- function(y, x, th = 10) {
  #bin.width <- var(x) * bw.SJ(x)
  #In some cases bw.SJ() failed with "sample too sparse to compute TD" due to gene(s) with infinite values
  gp <- names(x[which(!is.finite(x))])
  if(length(gp) == 0){
    bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
  }else{
    initg <- names(y)
    x <- x[is.finite(x)]
    y <- y[which(!names(y) %in% gp)]
    bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
  }
  eps <- .Machine$double.eps * 10
  breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
  breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
  score1 <- robust_scale_binned(y, x, breaks1)
  score2 <- robust_scale_binned(y, x, breaks2)
  res <- pmin(abs(score1), abs(score2)) > th
  if(length(gp)>0){
    #Add gp as outliers
    names(res) <- names(y)
    len <- length(res)
    gpres <- rep(TRUE,length(gp))
    names(gpres) <- gp
    res <- c(res,gpres)
    #Order the vector
    res <- res[order(match(names(res),initg))]
    #Remove the names
    res <- unname(res)
    res <- as.logical(res)
  }
  return(res)
}

#' Robust scale using median and mad per bin
#'
#' @param y Numeric vector
#' @param x Numeric vector
#' @param breaks Numeric vector of breaks
#'
#' @return Numeric vector of scaled score
#'
#' @importFrom stats aggregate
#'
robust_scale_binned <- function(y, x, breaks) {
  bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
  tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
  score <- rep(0, length(x))
  o <- order(bins)
  if (inherits(x = tmp$x, what = 'list')) {
    score[o] <- unlist(tmp$x)
  } else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score)
}

#' Robust scale using median and mad
#'
#' @param x Numeric
#'
#' @return Numeric
#'
#' @importFrom stats median mad
#'
robust_scale <- function(x) {
  return((x - median(x)) / (mad(x) + .Machine$double.eps))
}

pearson_residual <- function(y, mu, theta, min_var = -Inf) {
  model_var <- mu + mu^2 / theta
  model_var[model_var < min_var] <- min_var
  return((y - mu) / sqrt(model_var))
}

sq_deviance_residual <- function(y, mu, theta, wt=1) {
  2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
}

deviance_residual <- function(y, mu, theta, wt=1) {
  r <- 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
  sqrt(r) * sign(y - mu)
}

#' Return Pearson or deviance residuals of regularized models
#'
#' @param vst_out The output of a vst run
#' @param umi The UMI count matrix that will be used
#' @param residual_type What type of residuals to return; can be 'pearson' or 'deviance'; default is 'pearson'
#' @param res_clip_range Numeric of length two specifying the min and max values the results will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is vst_out$arguments$min_variance
#' @param cell_attr Data frame of cell meta data
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A matrix of residuals
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' pearson_res <- get_residuals(vst_out, pbmc)
#' deviance_res <- get_residuals(vst_out, pbmc, residual_type = 'deviance')
#' }
#'
get_residuals <- function(vst_out, umi, residual_type = 'pearson',
                          res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                          min_variance = vst_out$arguments$min_variance,
                          cell_attr = vst_out$cell_attr, bin_size = 256,
                          verbosity = vst_out$arguments$verbosity, 
                          verbose = NULL, show_progress = NULL) {
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }
  
  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  if (verbosity > 0) {
    message('Calculating residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res[genes_bin, ] <- switch(residual_type,
                               'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_variance),
                               'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                               stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment')
    )
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  res[res < res_clip_range[1]] <- res_clip_range[1]
  res[res > res_clip_range[2]] <- res_clip_range[2]
  return(res)
}

#' Return variance of residuals of regularized models
#'
#' This never creates the full residual matrix and can be used to determine highly variable genes.
#'
#' @param vst_out The output of a vst run
#' @param umi The UMI count matrix that will be used
#' @param residual_type What type of residuals to return; can be 'pearson' or 'deviance'; default is 'pearson'
#' @param res_clip_range Numeric of length two specifying the min and max values the residuals will be clipped to; default is c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
#' @param min_variance Lower bound for the estimated variance for any gene in any cell when calculating pearson residual; default is vst_out$arguments$min_variance
#' @param cell_attr Data frame of cell meta data
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A vector of residual variances (after clipping)
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' res_var <- get_residual_var(vst_out, pbmc)
#' }
#'
get_residual_var <- function(vst_out, umi, residual_type = 'pearson',
                             res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                             min_variance = vst_out$arguments$min_variance,
                             cell_attr = vst_out$cell_attr, bin_size = 256,
                             verbosity = vst_out$arguments$verbosity, 
                             verbose = NULL, show_progress = NULL) {
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }
  
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  model_pars <- vst_out$model_pars_fit
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }
  
  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars)]
  if (verbosity > 0) {
    message('Calculating variance for residuals of type ', residual_type, ' for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    y <- as.matrix(umi[genes_bin, , drop=FALSE])
    res_mat <- switch(residual_type,
                      'pearson' = pearson_residual(y, mu, model_pars[genes_bin, 'theta'], min_var = min_variance),
                      'deviance' = deviance_residual(y, mu, model_pars[genes_bin, 'theta']),
                      stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment'))
    res_mat[res_mat < res_clip_range[1]] <- res_clip_range[1]
    res_mat[res_mat > res_clip_range[2]] <- res_clip_range[2]
    res[genes_bin] <- row_var(res_mat)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  return(res)
}

#' Return average variance under negative binomial model
#'
#' This is based on the formula var = mu + mu^2 / theta
#'
#' @param vst_out The output of a vst run
#' @param cell_attr Data frame of cell meta data
#' @param use_nonreg Use the non-regularized parameter estimates; boolean; default is FALSE
#' @param bin_size Number of genes to put in each bin (to show progress)
#' @param verbosity An integer specifying whether to show only messages (1), messages and progress bars (2) or nothing (0) while the function is running; default is 2
#' @param verbose Deprecated; use verbosity instead
#' @param show_progress Deprecated; use verbosity instead
#'
#' @return A named vector of variances (the average across all cells), one entry per gene.
#'
#' @examples
#' \donttest{
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' res_var <- get_model_var(vst_out)
#' }
#'
get_model_var <- function(vst_out, cell_attr = vst_out$cell_attr, use_nonreg = FALSE,
                          bin_size = 256, verbosity = 2,
                          verbose = NULL, show_progress = NULL) {
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }
  
  regressor_data <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str)), cell_attr)
  if (use_nonreg) {
    model_pars <- vst_out$model_pars
  } else {
    model_pars <- vst_out$model_pars_fit
  }
  if (!is.null(dim(vst_out$model_pars_nonreg))) {
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', vst_out$model_str_nonreg)), cell_attr)
    regressor_data <- cbind(regressor_data, regressor_data_nonreg)
    model_pars <- cbind(vst_out$model_pars_fit, vst_out$model_pars_nonreg)
  }
  
  genes <- rownames(model_pars)
  if (verbosity > 0) {
    message('Calculating model variance for ', length(genes), ' genes')
  }
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 1) {
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  res <- matrix(NA_real_, length(genes))
  names(res) <- genes
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    mu <- exp(tcrossprod(model_pars[genes_bin, -1, drop=FALSE], regressor_data))
    model_var = mu + mu^2 / model_pars[genes_bin, 'theta']
    res[genes_bin] <- rowMeans(model_var)
    if (verbosity > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbosity > 1) {
    close(pb)
  }
  return(res)
}


############

vst.query <- function(umi,
                      cell_attr = NULL,
                      latent_var = c('log_umi'),
                      batch_var = NULL,
                      latent_var_nonreg = NULL,
                      n_genes = 2000,
                      n_cells = NULL,
                      method = 'poisson',
                      do_regularize = TRUE,
                      theta_regularization = 'od_factor',
                      res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
                      bin_size = 500,
                      min_cells = 5,
                      residual_type = 'pearson',
                      return_cell_attr = FALSE,
                      return_gene_attr = TRUE,
                      return_corrected_umi = FALSE,
                      min_variance = -Inf,
                      bw_adjust = 3,
                      gmean_eps = 1,
                      theta_estimation_fun = 'theta.ml',
                      theta_given = NULL,
                      verbosity = 2,
                      verbose = NULL,
                      show_progress = NULL,bg.model.pars) {
  arguments <- as.list(environment())
  arguments <- arguments[!names(arguments) %in% c("umi", "cell_attr")]
  
  # Take care of deprecated arguments
  if (!is.null(verbose)) {
    warning("The 'verbose' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    verbosity <- as.numeric(verbose)
  }
  if (!is.null(show_progress)) {
    warning("The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)", immediate. = TRUE, call. = FALSE)
    if (show_progress) {
      verbosity <- 2
    } else {
      verbosity <- min(verbosity, 1)
    }
  }
  
  # Check for suggested package
  if (method == "glmGamPoi") {
    glmGamPoi_check <- requireNamespace("glmGamPoi", quietly = TRUE)
    if (!glmGamPoi_check){
      stop('Please install the glmGamPoi package. See https://github.com/const-ae/glmGamPoi for details.')
    }
  }
  
  # Special case offset model - override most parameters
  if (startsWith(x = method, prefix = 'offset')) {
    cell_attr <- NULL
    latent_var <- c('log_umi')
    batch_var <- NULL
    latent_var_nonreg <- NULL
    n_genes <- NULL
    n_cells <- NULL
    do_regularize <- FALSE
    if (is.null(theta_given)) {
      theta_given <- 100
    } else {
      theta_given <- theta_given[1]
    }
  }
  
  
  times <- list(start_time = Sys.time())
  
  cell_attr <- make_cell_attr(umi, cell_attr, latent_var, batch_var, latent_var_nonreg, verbosity)
  if (!is.null(batch_var)) {
    cell_attr[, batch_var] <- as.factor(cell_attr[, batch_var])
    batch_levels <- levels(cell_attr[, batch_var])
  }
  
  # we will generate output for all genes detected in at least min_cells cells
  # but for the first step of parameter estimation we might use only a subset of genes
  genes_cell_count <- rowSums(umi >= 0.01)
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))
  
  if (!do_regularize && !is.null(n_genes)) {
    if (verbosity > 0) {
      message('do_regularize is set to FALSE, will use all genes')
    }
    n_genes <- NULL
  }
  
  if (!is.null(n_cells) && n_cells < ncol(umi)) {
    # downsample cells to speed up the first step
    cells_step1 <- sample(x = colnames(umi), size = n_cells)
    if (!is.null(batch_var)) {
      dropped_batch_levels <- setdiff(batch_levels, levels(droplevels(cell_attr[cells_step1, batch_var])))
      if (length(dropped_batch_levels) > 0) {
        stop('Dropped batch levels ', dropped_batch_levels, ', set n_cells higher')
      }
    }
    genes_cell_count_step1 <- rowSums(umi[, cells_step1] > 0)
    genes_step1 <- rownames(umi)[genes_cell_count_step1 >= min_cells]
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  } else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean
  }
  
  data_step1 <- cell_attr[cells_step1, , drop = FALSE]
  
  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    # density-sample genes to speed up the first step
    log_gmean_dens <- density(x = genes_log_gmean_step1, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = log_gmean_dens$x, y = log_gmean_dens$y, xout = genes_log_gmean_step1)$y + .Machine$double.eps)
    genes_step1 <- sample(x = genes_step1, size = n_genes, prob = sampling_prob)
    genes_log_gmean_step1 <- log10(row_gmean(umi[genes_step1, cells_step1], eps = gmean_eps))
  }
  
  if (!is.null(batch_var)) {
    model_str <- paste0('y ~ (', paste(latent_var, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
  } else {
    model_str <- paste0('y ~ ', paste(latent_var, collapse = ' + '))
  }
  
  bin_ind <- ceiling(x = 1:length(x = genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  if (verbosity > 0) {
    message('Variance stabilizing transformation of count matrix of size ', nrow(umi), ' by ', ncol(umi))
    message('Model formula is ', model_str)
  }
  
  times$get_model_pars = Sys.time()
  model_pars <- bg.model.pars
  # make sure theta is not too small
  min_theta <- 1e-7
  if (any(model_pars[, 'theta'] < min_theta)) {
    if (verbosity > 0) {
      msg <- sprintf('There are %d estimated thetas smaller than %g - will be set to %g', sum(model_pars[, 'theta'] < min_theta), min_theta, min_theta)
      message(msg)
    }
    model_pars[, 'theta'] <- pmax(model_pars[, 'theta'], min_theta)
  }
  
  model_pars <- model_pars[names(genes_log_gmean_step1),]
  
  times$reg_model_pars = Sys.time()
  if (do_regularize) {
    model_pars_fit <- reg_model_pars(model_pars = model_pars, genes_log_gmean_step1 = genes_log_gmean_step1,
                                     genes_log_gmean = genes_log_gmean,cell_attr =  cell_attr,
                                     batch_var = batch_var,cells_step1 =  cells_step1,
                                     genes_step1 =  genes_step1,umi = umi,bw_adjust =  bw_adjust,gmean_eps =  gmean_eps,
                                     theta_regularization = theta_regularization,verbosity = verbosity)
    model_pars_outliers <- attr(model_pars_fit, 'outliers')
  } else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }
  
  # use all fitted values in NB model
  regressor_data <- model.matrix(as.formula(gsub('^y', '', model_str)), cell_attr)
  
  if (!is.null(latent_var_nonreg)) {
    if (verbosity > 0) {
      message('Estimating parameters for following non-regularized variables: ', latent_var_nonreg)
    }
    if (!is.null(batch_var)) {
      model_str_nonreg <- paste0('y ~ (', paste(latent_var_nonreg, collapse = ' + '), ') : ', batch_var, ' + ', batch_var, ' + 0')
    } else {
      model_str_nonreg <- paste0('y ~ ', paste(latent_var_nonreg, collapse = ' + '))
    }
    
    times$get_model_pars_nonreg = Sys.time()
    model_pars_nonreg <- get_model_pars_nonreg(genes, bin_size, model_pars_fit, regressor_data, umi, model_str_nonreg, cell_attr, verbosity)
    
    regressor_data_nonreg <- model.matrix(as.formula(gsub('^y', '', model_str_nonreg)), cell_attr)
    model_pars_final <- cbind(model_pars_fit, model_pars_nonreg)
    regressor_data_final <- cbind(regressor_data, regressor_data_nonreg)
    #model_pars_final[, '(Intercept)'] <- model_pars_final[, '(Intercept)'] + model_pars_nonreg[, '(Intercept)']
    #model_pars_final <- cbind(model_pars_final, model_pars_nonreg[, -1, drop=FALSE])
    # model_str <- paste0(model_str, gsub('^y ~ 1', '', model_str2))
  } else {
    model_str_nonreg <- ''
    model_pars_nonreg <- c()
    model_pars_final <- model_pars_fit
    regressor_data_final <- regressor_data
  }
  
  times$get_residuals = Sys.time()
  if (!residual_type == 'none') {
    if (verbosity > 0) {
      message('Second step: Get residuals using fitted parameters for ', length(x = genes), ' genes')
    }
    bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
    max_bin <- max(bin_ind)
    if (verbosity > 1) {
      pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
    }
    res <- matrix(NA_real_, length(genes), nrow(regressor_data_final), dimnames = list(genes, rownames(regressor_data_final)))
    for (i in 1:max_bin) {
      genes_bin <- genes[bin_ind == i]
      mu <- exp(tcrossprod(model_pars_final[genes_bin, -1, drop=FALSE], regressor_data_final))
      y <- as.matrix(umi[genes_bin, , drop=FALSE])
      res[genes_bin, ] <- switch(residual_type,
                                 'pearson' = pearson_residual(y, mu, model_pars_final[genes_bin, 'theta'], min_var = min_variance),
                                 'deviance' = deviance_residual(y, mu, model_pars_final[genes_bin, 'theta']),
                                 stop('residual_type ', residual_type, ' unknown - only pearson and deviance supported at the moment')
      )
      if (verbosity > 1) {
        setTxtProgressBar(pb, i)
      }
    }
    if (verbosity > 1) {
      close(pb)
    }
  } else {
    if (verbosity > 0) {
      message('Skip calculation of full residual matrix')
    }
    res <- matrix(data = NA, nrow = 0, ncol = 0)
  }
  
  rv <- list(y = res,
             model_str = model_str,
             model_pars = model_pars,
             model_pars_outliers = model_pars_outliers,
             model_pars_fit = model_pars_fit,
             model_str_nonreg = model_str_nonreg,
             model_pars_nonreg = model_pars_nonreg,
             arguments = arguments,
             genes_log_gmean_step1 = genes_log_gmean_step1,
             cells_step1 = cells_step1,
             cell_attr = cell_attr)
  rm(res)
  gc(verbose = FALSE)
  
  times$correct_umi = Sys.time()
  if (return_corrected_umi) {
    if (residual_type != 'pearson') {
      message("Will not return corrected UMI because residual type is not set to 'pearson'")
    } else {
      rv$umi_corrected <- sctransform::correct(rv, do_round = TRUE, do_pos = TRUE,
                                               verbosity = verbosity)
      rv$umi_corrected <- as(object = rv$umi_corrected, Class = 'dgCMatrix')
    }
  }
  
  rv$y[rv$y < res_clip_range[1]] <- res_clip_range[1]
  rv$y[rv$y > res_clip_range[2]] <- res_clip_range[2]
  
  if (!return_cell_attr) {
    rv[['cell_attr']] <- NULL
  }
  
  times$get_gene_attr = Sys.time()
  if (return_gene_attr) {
    if (verbosity > 0) {
      message('Calculating gene attributes')
    }
    gene_attr <- data.frame(
      detection_rate = genes_cell_count[genes] / ncol(umi),
      gmean = 10 ^ genes_log_gmean,
      variance = row_var(umi))
    if (ncol(rv$y) > 0) {
      gene_attr$residual_mean = rowMeans(rv$y)
      gene_attr$residual_variance = row_var(rv$y)
    }
    # Special case offset model - also calculate arithmetic mean
    if (startsWith(x = method, prefix = 'offset')) {
      gene_attr$amean <- rowMeans(umi)
    }
    
    rv[['gene_attr']] <- gene_attr
  }
  
  if (verbosity > 0) {
    message('Wall clock passed: ', capture.output(print(Sys.time() - times$start_time)))
  }
  times$done = Sys.time()
  rv$times <- times
  return(rv)
}
