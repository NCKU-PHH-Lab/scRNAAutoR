MetaSummary = function(scRNA_SeuObj.list,  scRNA.SeuObj,
                       scRNA_SeuObj_QC.list, scRNA_Ori.SeuObj,
                       SavePath = "", projectName = ""){
  
  Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
  colnames(Meta.df) <- c("Folder","Cell_Num","Gene_Num")
  
  ## Before QC
  for (i in 1:length(scRNA_SeuObj.list)) {
    Meta.df[i,1] <- names(scRNA_SeuObj.list[i])
    Meta.df[i,2] <- ncol(scRNA_SeuObj.list[[i]]@assays[["RNA"]]@counts)
    Meta.df[i,3] <- nrow(scRNA_SeuObj.list[[i]]@assays[["RNA"]]@counts)
    
  }
  
  # Summary to Meta table
  Meta.df[i+1,1] <- c("Summary")
  Meta.df[i+1,2] <- ncol(scRNA_Ori.SeuObj@assays[["RNA"]]@counts)
  Meta.df[i+1,3] <- nrow(scRNA_Ori.SeuObj@assays[["RNA"]]@counts)
  
  ## After QC
  for (j in 1:length(scRNA_SeuObj_QC.list)) {
    Meta.df[i+j+1,1] <- paste0(names(scRNA_SeuObj_QC.list[j]),".QC")
    Meta.df[i+j+1,2] <- ncol(scRNA_SeuObj_QC.list[[j]]@assays[["RNA"]]@counts)
    Meta.df[i+j+1,3] <- nrow(scRNA_SeuObj_QC.list[[j]]@assays[["RNA"]]@counts)
    
  }
  
  # Summary to Meta table
  Meta.df[i+j+2,1] <- c("Summary.QC")
  Meta.df[i+j+2,2] <- ncol(scRNA.SeuObj@assays[["RNA"]]@counts)
  Meta.df[i+j+2,3] <- nrow(scRNA.SeuObj@assays[["RNA"]]@counts)
  
  
  rm(i,j)
  
  Meta.df <- left_join(Meta.df,list_files.df, by = "Folder")
  Meta.df[is.na(Meta.df)] <- ""
  
  ## Export data        
  write.table( Meta.df ,
               file = paste0(SavePath,projectName,"_CellCount_Meta.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )
  
  
  return(Meta.df)
}