ReadscRNA = function( InputFolder, list_files.df, Mode="10x" ,ProjectName="") # Mode=c("10x","Exp")
  {

    if(Mode=="10x"){
      ## Read 10x files
      scRNA_SeuObj.list <- list()
      for(i in 1:nrow(list_files.df)){
        Folder <- list_files.df$Folder[i]
        Data.dgCMatrix <- Read10X(data.dir = paste0(InputFolder,"/", Folder, "/monocle/outs/filtered_gene_bc_matrices/mm10"))
        Data.SeuObj <- CreateSeuratObject(counts = Data.dgCMatrix, 
                                          project = ProjectName, 
                                          min.cells = 3, min.features = 200)
        
        for (j in 1:(ncol(list_files.df)-1)) {
          Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1], 
                                                                       times=length(Data.SeuObj@meta.data[["orig.ident"]]))
        }
        
        scRNA_SeuObj.list[[i]] <- Data.SeuObj
        names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]
      }
      rm(i,j,Folder,Data.dgCMatrix,Data.SeuObj)
      
    }else{
      
      ## Expression matrix
      ## GSE103322 HNSC
      scRNA_SeuObj.list <- list()
      for(i in 1:nrow(list_files.df)){
        Folder <- list_files.df$Folder[i]
        GeneExp.df <- read.table(paste0(InputFolder,"/", Folder, "/Exp.tsv"),
                                     header=T, row.names = 1, sep="\t")
        Data.SeuObj <- CreateSeuratObject(counts = GeneExp.df,
                                          project = ProjectName,
                                          min.cells = 3, min.features = 200)
        
        for (j in 1:(ncol(list_files.df)-1)) {
          Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1], 
                                                                       times=length(Data.SeuObj@meta.data[["orig.ident"]]))
        }
        
        scRNA_SeuObj.list[[i]] <- Data.SeuObj
        names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]
        
        scAnno.df <- read.table(paste0(InputFolder,"/", Folder, "/Anno.tsv"),
                                     header=T, row.names = 1, sep="\t")
        # scRNA.SeuObj@meta.data[["sample"]] <- GeneExp.df[5,] %>% as.character()
        
        for (k in 1:nrow(scAnno.df)) {
          scRNA_SeuObj.list[[i]]@meta.data[[row.names(scAnno.df)[k]]] <- c(scAnno.df[k,] %>% as.character())
        }
        
      }
      rm(i,j,Folder,GeneExp.df,Data.SeuObj)
 
    }
    return(scRNA_SeuObj.list)
  }
    