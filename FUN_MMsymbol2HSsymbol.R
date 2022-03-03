##### Function design reference#####
# Reference_Alias2HSsymbol: https://www.biostars.org/p/14971/
# Reference_Alias2MMsymbol: http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/alias2Symbol.html
# Reference_HSsymbol2MMsymbol: https://www.biostars.org/p/149115/


##### Alias to Human gene symbol (HGNC) #####
Alias2HSsymbol <- function(df,alias_col){
  
  require(org.Hs.eg.db)
  
  # create new column
  df <- data.frame(df,"HS-symbol"=rep(NA,nrow(df)))
  cat("Create new column: OK","\n")
  
  # crate alias set
  query <- toupper(df[[alias_col]])
  cat("Create alias query: OK","\n")
  
  # use sql to get alias table and gene_info table (contains the symbols)
  
  ## first open the database connection
  dbCon <- org.Hs.eg_dbconn()
  ## write your SQL query
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  ## execute the query on the database
  aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
  cat("Data query: OK","\n")
  
  # transformation
  cat("Start alias transformation","\n")
  count <- 1
  for (i in 1:length(query)){
    Hmsymbol <- aliasSymbol[which(aliasSymbol[,2] == query[i]),5]
    n_Hmsymbol <- length(Hmsymbol)
    cat(" Alias",i,"\n")
    
    if(n_Hmsymbol == 1){
      df[count,ncol(df)] <- Hmsymbol
      count <- count+1
      cat("  Symbol","1","\n")
    }else if(n_Hmsymbol > 1){
      ## Expend the row that having multiple symbol
      working.df <- df[c(1:count-1),]
      for (j in 1:n_Hmsymbol){
        working.df <- rbind(working.df,df[count,])
      }
      working.df <- rbind(working.df,df[-c(1:count),])
      df <- working.df
      rm(working.df)
      ## Insert the gene symbol
      for (k in 1:n_Hmsymbol) {
        df[count+k-1,ncol(df)] <- Hmsymbol[k]
      }
      count <- count+k
      cat("  Symbol",k,"\n")
    }else{
      df[count,ncol(df)] <- 0
      count <- count+1
      cat("  Symbol not found","\n")
    }
  }
  row.names(df) <- c(1:nrow(df))
  cat("Finish")
  return(df)
}

##### Alias to Mouse gene symbol (MGI) #####
Alias2MMsymbol <- function(df,alias_col){
  
  require(limma)
  
  # create new column
  df <- data.frame(df, "MM-symbol"=rep(NA,nrow(df)))
  cat("Create new column: OK","\n")
  
  # crate alias set
  query <- df[[alias_col]]
  cat("Create alias query: OK","\n")
  
  # transformation
  cat("Start alias transformation","\n")
  count <- 1
  for (i in 1:length(query)) {
    mmsymbol <- alias2Symbol(query[i], species = "Mm")
    n_mmsymbol <- length(mmsymbol)
    cat(" Alias",i,"\n")
    
    if(n_mmsymbol==1){
      df[count,ncol(df)] <- mmsymbol
      count <- count+1
      cat("  Symbol","1","\n")
    }else if(n_mmsymbol>1){
      working.df <- df[c(1:count-1),]
      
      ## Expend the row that having multiple symbol
      for(j in 1: n_mmsymbol){
        working.df <- rbind(working.df,df[count,])
      }
      working.df <- rbind(working.df,df[-c(1:count),])
      df <- working.df
      rm(working.df)
      
      ## Insert the gene symbol
      for (k in 1:n_mmsymbol) {
        df[count+k-1,ncol(df)] <- mmsymbol[k]
      }
      count <- count+k
      cat("  Symbol",k,"\n")
    }else{
      df[count,ncol(df)] <- 0
      count <- count+1
      cat("  Symbol not found","\n")
    }  
  }
  row.names(df) <- c(1:nrow(df))
  cat("Finish")
  return(df)
}

##### Human gene symbol(HGNC) to Mouse gene symbol (MGI) #####
HSsymbol2MMsymbol<- function(df,HS_col){
  
  require(epanetReader)
  
  # download the dataset
  if(!file.exists("HOM_MouseHumanSequence.rpt")){
    download.file("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", destfile = "HOM_MouseHumanSequence.rpt")
  }
  MouseHumanSequence.df <- read.csv(file=paste0(getwd(),"/HOM_MouseHumanSequence.rpt"), sep = "\t")
  
  # create new column
  df <- data.frame(df, "MM-symbol"=rep(NA,nrow(df)))
  cat("Create new column: OK","\n")
  
  # crate alias set
  query <- df[[HS_col]]
  cat("Create alias query: OK","\n")
  
  # transformation
  cat("Start alias transformation","\n")
  count <- 1
  for (i in 1:length(query)) {
    
    HSID <- as.numeric(row.names(MouseHumanSequence.df[MouseHumanSequence.df$Symbol==query[i],] ))
    MMID <- c()
    if(length(HSID)>0){
      for (l in 1:length(HSID)) {
        n <- HSID[l]
        while(MouseHumanSequence.df$Common.Organism.Name[n]!="mouse, laboratory"){
          n <- n-1
        }
        MMID[l] <- n
      }
      mmsymbol <- MouseHumanSequence.df$Symbol[MMID]
    }else{
      mmsymbol <- c()
    }
    
    n_mmsymbol <- length(mmsymbol)
    
    cat(" HS",i,"\n")
    if(n_mmsymbol==1){
      df[count,ncol(df)] <- mmsymbol
      count <- count+1
      cat("  MM","1","\n")
    }else if(n_mmsymbol>1){
      working.df <- df[c(1:count-1),]
      
      ## Expend the row that having multiple symbol
      for(j in 1: n_mmsymbol){
        working.df <- rbind(working.df,df[count,])
      }
      working.df <- rbind(working.df,df[-c(1:count),])
      df <- working.df
      rm(working.df)
      
      ## Insert the gene symbol
      for (k in 1:n_mmsymbol) {
        df[count+k-1,ncol(df)] <- mmsymbol[k]
      }
      count <- count+k
      cat("  MM",k,"\n")
    }else{
      df[count,ncol(df)] <- 0
      count <- count+1
      cat("  MM not found","\n")
    }  
  }
  row.names(df) <- c(1:nrow(df))
  cat("Finish")
  return(df)
}