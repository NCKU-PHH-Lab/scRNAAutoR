AnnoSummary = function(scRNA.SeuObj,  list_files.df, scAnno.df,
                       ClassSet=1, ClassSet2 = ClassSet2){

  ## Annotation Summary Table
    for (i in 1:(ncol(list_files.df)+nrow(scAnno.df)-1)) {
      if(i==1){
        Anno.df <- data.frame(scRNA.SeuObj@meta.data[[3+i]])
        colnames(Anno.df)[i] <- colnames(list_files.df)[i+1]
      }else{
        Anno_S.df <- data.frame(scRNA.SeuObj@meta.data[[3+i]])
        Anno.df <- data.frame(Anno.df, Anno_S.df)
        colnames(Anno.df)[i] <- colnames(list_files.df)[i+1]
      }
    }
    rm(i,Anno_S.df)
    
    Anno.df <- data.frame(Anno.df, celltype = scRNA.SeuObj@meta.data[["celltype"]])
  
  ## Annotation Frequency
    # All
    Anno_Freq_df.lt <- list()
    
    for (i in 1:(ncol(Anno.df))) {
      Anno_Freq_df.lt[[i]] <- table(Anno.df[,i]) %>% as.data.frame()
      names(Anno_Freq_df.lt)[[i]] <- paste0("Freq_",colnames(Anno.df)[i])
      
    }
    rm(i)
    
    
    # Group by sample
    # assign(colnames(Anno.df)[1], Anno.df[,1] %>% unique()) 
    Anno_Tar.set <- Anno.df[, ClassSet] %>% unique()
    Anno_Tar.Num <- Anno_Tar.set %>% length()
    Anno_Tar_df.lt <- list()
    
    for (i in 1:Anno_Tar.Num) {
      Anno_Tar_df.lt[[i]] <- Anno.df[Anno.df[, ClassSet]==Anno_Tar.set[i],]
      names(Anno_Tar_df.lt)[[i]] <-  paste0("State_",Anno_Tar.set[i])
    }
    rm(i)
    
    Anno_Freq_Tar_df.lt <- list()
    
    for (j in 1:length(Anno_Tar_df.lt)) {
      Anno_Freq_Tar_df.lt[[j]] <- list()
      names(Anno_Freq_Tar_df.lt)[[j]] <- paste0(names(Anno_Tar_df.lt)[j])
      
      for (i in 1:(ncol(Anno_Tar_df.lt[[j]]))) {
        Anno_Freq_Tar_df.lt[[j]][[i]] <- table(Anno_Tar_df.lt[[j]][,i]) %>% as.data.frame()
        names(Anno_Freq_Tar_df.lt[[j]][[i]])[1] <- colnames(Anno.df)[i]
        Anno_Freq_Tar_df.lt[[j]][[i]] <- data.frame(Type=paste0(Anno_Tar_df.lt[[j]][1,1]),
                                                    Anno_Freq_Tar_df.lt[[j]][[i]])
        colnames(Anno_Freq_Tar_df.lt[[j]][[i]])[1] <- colnames(list_files.df)[ClassSet+1]
        Anno_Freq_Tar_df.lt[[j]][[i]]$Percent <- Anno_Freq_Tar_df.lt[[j]][[i]]$Freq/sum(Anno_Freq_Tar_df.lt[[j]][[i]]$Freq)
        
        names(Anno_Freq_Tar_df.lt[[j]])[[i]] <- paste0("Freq_",colnames(Anno_Tar_df.lt[[j]])[[i]])
        
      }
      if(j==1){
        Freq_All.df <- Anno_Freq_Tar_df.lt[[j]][[i]]
      }else{
        Freq_All.df <- rbind(Freq_All.df, Anno_Freq_Tar_df.lt[[j]][[i]])
      }
      
    }
    rm(i,j)
  
  # Combind all count of sample
    Freq_All.df <- data.frame(Index = row.names(Freq_All.df),Freq_All.df )
    colnames(Freq_All.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")

  #### All type compare to Combine Sex ####  
    ##
    Freq_All_Cla.lt <- list() 
    
    ClassSet2.set <- Anno.df[,ClassSet2] %>% unique()
    Freq_All_Cla.lt[[paste0("Anno_Cla",1)]] <- Anno.df[Anno.df[,ClassSet2] == ClassSet2.set[1],]
    Freq_All_Cla.lt[[paste0("Anno_Cla",2)]] <- Anno.df[Anno.df[,ClassSet2] == ClassSet2.set[2],]
    
    
    # Count EO_CT
    Freq_All_Cla.lt[[paste0("Freq_Cla",1)]] <- table(Freq_All_Cla.lt[[paste0("Anno_Cla",1)]]$celltype) %>% as.data.frame()
    Freq_All_Cla.lt[[paste0("Freq_Cla",1)]] <- data.frame(Type=ClassSet2.set[1],Freq_All_Cla.lt[[paste0("Freq_Cla",1)]])
    Freq_All_Cla.lt[[paste0("Freq_Cla",1)]]$Percent <- Freq_All_Cla.lt[[paste0("Freq_Cla",1)]]$Freq/sum(Freq_All_Cla.lt[[paste0("Freq_Cla",1)]]$Freq)
    colnames(Freq_All_Cla.lt[[paste0("Freq_Cla",1)]]) <- c("Sample","celltype","Freq","Percent")
    
    # Count LO_CT
    Freq_All_Cla.lt[[paste0("Freq_Cla",2)]] <- table(Freq_All_Cla.lt[[paste0("Anno_Cla",2)]]$celltype) %>% as.data.frame()
    Freq_All_Cla.lt[[paste0("Freq_Cla",2)]] <- data.frame(Type=ClassSet2.set[2],Freq_All_Cla.lt[[paste0("Freq_Cla",2)]])
    Freq_All_Cla.lt[[paste0("Freq_Cla",2)]]$Percent <- Freq_All_Cla.lt[[paste0("Freq_Cla",2)]]$Freq/sum(Freq_All_Cla.lt[[paste0("Freq_Cla",2)]]$Freq)
    colnames(Freq_All_Cla.lt[[paste0("Freq_Cla",2)]]) <- c("Sample","celltype","Freq","Percent")
    
    # Combind all count of sample
    Freq_All_Cla.lt[["Freq_All_Cla.df"]]  <- rbind(Anno_Freq_Tar_df.lt[[1]][["Freq_celltype"]],
                                                   Anno_Freq_Tar_df.lt[[2]][["Freq_celltype"]],
                                                   Anno_Freq_Tar_df.lt[[3]][["Freq_celltype"]],
                                                   Anno_Freq_Tar_df.lt[[4]][["Freq_celltype"]],
                                                   Freq_All_Cla.lt[[paste0("Freq_Cla",1)]],
                                                   Freq_All_Cla.lt[[paste0("Freq_Cla",2)]])
    
  
  
    Freq_All_Cla.lt[["Freq_All_Cla.df"]] <- data.frame(Index = row.names(Freq_All_Cla.lt[["Freq_All_Cla.df"]]),Freq_All_Cla.lt[["Freq_All_Cla.df"]] )
    colnames(Freq_All_Cla.lt[["Freq_All_Cla.df"]]) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")
    
    
    
    Freq_All_Cla.lt[["Freq_All_Cla.df"]]$Percent <- as.numeric(Freq_All_Cla.lt[["Freq_All_Cla.df"]]$Percent)
  
  
  
#### Export data #####  
  AnnoSummary.lt <- list()
  AnnoSummary.lt[["Anno.df"]] <- Anno.df
  AnnoSummary.lt[["Anno_Tar.set"]] <- Anno_Tar.set
  AnnoSummary.lt[["Anno_Tar.Num"]] <- Anno_Tar.Num
  AnnoSummary.lt[["Anno_Tar_df.lt"]] <- Anno_Tar_df.lt
  AnnoSummary.lt[["Anno_Freq_Tar_df.lt"]] <- Anno_Freq_Tar_df.lt
  AnnoSummary.lt[["Anno_Tar_df.lt"]] <- Anno_Tar_df.lt
  AnnoSummary.lt[["Freq_All.df"]] <- Freq_All.df
  AnnoSummary.lt[["Freq_All_Cla.lt"]] <- Freq_All_Cla.lt
  
  return(AnnoSummary.lt)
}  
  