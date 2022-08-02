#' scRefBool: quantifies different levels of gene expression including low, medium and high
#'
#' @param data single cell expression matrix
#' @param metadata corresponding metadata (two columns: cell names, cell type)
#' @param org organism (human or mouse)
#' @param dir directory used to save the results
#' @param file.name prefix for output files and folder
#' @param precBack precompiled background path to ScRbRef
#' @param discretize gene quantification in three levels of cell types
#' @param sig.frames construction of significant genes data frame 
#' @param ncores number of cores to use to run the code
#' @param fixseed seed used to discretize cell types (reproducbility)
#'
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import reshape2
#' @import tictoc
#' @import SingleCellExperiment
#' @import Matrix
#' @import scater
#' @import stringr
#' @import qlcMatrix
#' @import sctransform
## @import org.Mm.eg.db
## @import org.Hs.eg.db
#' @import scuttle
#' @import LaplacesDemon
#'
#' @export
#' 
#' @return NULL, all results are accessible in the dir specified
#' @author Céline Barlier, Kartikeya Singh, Sascha Jung
run_scRefBool <- function(data,metadata,org="human",dir,file.name,precBack=NA,discretize=TRUE,sig.frames=TRUE,ncores=detectCores()-2,fixseed=1234){

  # Changing column names of metadata for consistency
  if(length(colnames(metadata)) == 3){
    colnames(metadata) <- c("cell.id","tissue","cell.type")
  }else{
    colnames(metadata) <- c("cell.id","cell.type")
  }
  rownames(metadata) <- metadata$cell.id

  # Transforming matrix to sparse matrix if needed
  if(class(data)[1] == "matrix"){
    data <- as(data, "dgCMatrix")
  }
  
  # Identifying common cells between gene expression data frame and metadata - cleaning data and metadata
  common.cells <- intersect(metadata$cell.id,dimnames(data)[[2]])
  metadata <- metadata[which(metadata$cell.id %in% common.cells),]
  data <- data[,dimnames(data)[[2]] %in% common.cells]

  ################### BACKGROUND CONSTRUCTION/SELECTION #######################

  #If not pre-compiled background specified, build background
  if(is.na(precBack)){

    # Calculating number of cells of different cell types in data
    cell.type.freq <- as.data.frame(table(metadata$cell.type))
    # Need at least 50 cells by cell type
    cell.type.to.keep <- as.character(cell.type.freq[which(cell.type.freq$Freq >= 50),"Var1"])
    metadataBack <- metadata[which(metadata$cell.type %in% cell.type.to.keep),]
    dataBack <- data[,as.character(metadataBack$cell.id)]

    if(length(cell.type.to.keep) < 2){

      cat("Data contains less than 2 celltypes with at least 50 cells. Background reference cannot be created.")

    }else{

      removed.cell.type <- setdiff(cell.type.freq$Var1,cell.type.to.keep)

      #Create directory for results
      if(!dir.exists(paths = paste0(dir,"/Ref_",file.name))){
        dir.create(paste0(dir,"/Ref_",file.name))
      }
      
      #Create log file
      log.file.name <- paste0(dir,"/Ref_",file.name,"/Log_file_",paste(unlist(strsplit(x = as.character(Sys.time()),split = " ",fixed = 2)),collapse = "_"),".txt")
      write(x = paste0("Cell types used for Reference background :",paste(cell.type.to.keep,collapse = ",")),file = log.file.name)
      write(x = paste0("Cell types removed :",paste(removed.cell.type,collapse = ",")),file = log.file.name,append = T)

      setwd(dir)
      #Cleaning & Normalizing background data
      if(!file.exists(paste0("./Ref_",file.name,"/",file.name,"_NormData.rds"))){
        clean.data <- suppressWarnings(cleanData(data_counts = dataBack,dir.name=paste(dir,paste("/Ref_",file.name,sep=""),sep=""),file.id = file.name))
        norm.data <- suppressWarnings(scrb_norm(data_counts = clean.data,file.id = file.name,saveData = T,return = T))
      }else{
        norm.data <- readRDS(paste0("./Ref_",file.name,"/",file.name,"_NormData.rds"))
      }
      
      #Sampling, scaling & getting background thresholds
      if(!file.exists(paste0("./Ref_",file.name,"/",file.name,"_th_dist.rds"))){
        thrs.data <- GetBackground(data_norm=norm.data,dir_name=paste0(dir,"/Ref_",file.name),file.id=file.name,metadata=metadataBack,sample_n = 100,parallelize=T,ncores=ncores,num.boot.sample = 1000,fixseed = fixseed)
      }else{
        thrs.data <- readRDS(paste0("./Ref_",file.name,"/",file.name,"_th_dist.rds"))
      }
      
      setwd(paste0("./Ref_",file.name))
      #Normalization parameters used for the background (the algorithm will use the same to normalize the queries)
      bg.model.pars <- readRDS(file = paste0(file.name,"_bg_model_pars.rds"))
    }

  }else{

    #Select pre-compiled background
    thrs.data <- readRDS(list.files(path = precBack,pattern = "_th_dist.rds",full.names = T))
    bg.model.pars <- readRDS(list.files(path = precBack,pattern = "_bg_model_pars.rds",full.names = T))
    sfback <- read.delim(list.files(precBack,pattern = "sf",full.names = T),sep=" ",row.names = 1,header = F)
    bimod.genes <- readRDS(list.files(path = precBack,pattern = "_bimod.rds",full.names = T))
  }

  #############################################################################

  ################### GENE QUANTIFICATION ####################

  # Keep cell type/subtype query with at least 10 cells
  query.cell.type.freq <- as.data.frame(table(metadata$cell.type))
  query.cell.type.to.keep <- as.character(query.cell.type.freq[which(query.cell.type.freq$Freq >= 10),"Var1"])
  metadataQuery <- metadata[which(metadata$cell.type %in% query.cell.type.to.keep),]
  dataQuery <- data[,as.character(metadataQuery$cell.id)]
  
  if(length(query.cell.type.to.keep) == 0){
    cat("Data contains no cell (sub)type with at least 10 cells")
  }else{
    
    #Discretize cell type queries gene expression
    if(discretize){
      #Fix seed
      set.seed(fixseed)
      cat("Discretising data ...")
      lapply(X = query.cell.type.to.keep,FUN = function(pop){
        if(!file.exists(paste0("./Query_",pop,"/",pop,"_DiscretisedData.rds"))){
          cat("Celltype : ",pop,"\n")
          pop.meta <- metadataQuery[which(metadataQuery$cell.type == pop),]
          pop.data <- dataQuery[,as.character(pop.meta$cell.id)]
          #If no pre-compiled background is used
          if(is.na(precBack)){
            #Ensure background exists
            if(file.exists(paste0(dir,"/Ref_",file.name,"/",file.name,"_th_dist.rds"))){
              pop.scaled <- suppressWarnings(PrepQueryData(query_dat = pop.data,Ref_data = paste0(dir,"/Ref_",file.name),dir=paste0(dir,"/Ref_",file.name,"/Query_",pop),ncores=ncores,file.id = pop,saveData = T,return.data = T,bg.model.pars = bg.model.pars))
              pop.discrete <- Discretise_data(query_ScaledData = pop.scaled,file.id = pop,Ref_data = "./",saveData = T,return.data = F)
            }else{
              cat("No background reference found")
            }
          }else{
            setwd(dir)
            pop.scaled <- suppressWarnings(PrepQueryData(query_dat = pop.data,Ref_data = precBack,dir=paste0(dir,"/Query_",pop),sfback=sfback,ncores=ncores,file.id = pop,saveData = T,return.data = T,bg.model.pars = bg.model.pars))
            pop.discrete <- Discretise_data(query_ScaledData = pop.scaled,file.id = pop,Ref_data = precBack,thsback=thrs.data,saveData = T,return.data = F)
          }
        }
      })
    }
    
    #Significant levels of gene expression
    if(sig.frames){
       if(!is.na(precBack)){
          setwd(dir)
       }
      #If cell type discretized
      if(file.exists(paste0("./Query_",query.cell.type.to.keep[1],"/",query.cell.type.to.keep[1],"_DiscretisedData.rds"))){
        #Gene level by population
        dir.create("sig_frames")
        pop.sig.frames <- lapply(X = query.cell.type.to.keep,FUN = function(pop){
          cat("Cell type : ",pop,"\n")
          discrete.data <- readRDS(file = paste0("./Query_",pop,"/",pop,"_DiscretisedData.rds"))
          discrete.mat <- discrete.data$discretised_data
          rm(discrete.data)
          invisible(gc())
          sig.frame <- Create.sig.frame(discrete.mat = discrete.mat)
          sig.frame.expr.level <- discrete.expr.levels(tbl = sig.frame[[1]], zs = sig.frame[[2]])
          #Save
          write.table(x = sig.frame.expr.level,file = paste0("./sig_frames/",pop,"_sig_tbl.txt"),sep = "\t",row.names = F,quote = F)
        })
      }
    }
  }

  ##############################################################

  return(NULL)
}
