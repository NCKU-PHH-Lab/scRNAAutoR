## Function of count gene number of each cell type
CountGeneNum <- function(CCMar_SexCom_Split_Pos.list, mode=1) {

  if(mode==1){ # mode1 Count genes by different sub type
    for (i in c(1:length(CCMar_SexCom_Split_Pos.list))){
      CCMar_SexCom_Split_Pos_CT.list <- CCMar_SexCom_Split_Pos.list[i]
      names(CCMar_SexCom_Split_Pos_CT.list) <- names(CCMar_SexCom_Split_Pos.list[i])
      if(i == 1){
        CCMar_SexCom_Split_Pos_CT.df <- CountGene(CCMar_SexCom_Split_Pos_CT.list)
        colnames(CCMar_SexCom_Split_Pos_CT.df) <- c("Genes",names(CCMar_SexCom_Split_Pos.list[i]))
      }else{
        attach.df <- CountGene(CCMar_SexCom_Split_Pos_CT.list)
        if(ncol(attach.df) == 1 ){
          attach.df <-  as.data.frame(matrix(nrow=0,ncol=2)) 

          #attach.df = data.frame(Genes="",CellType="")
        }else{ }
        colnames(attach.df) <- c("Genes",names(CCMar_SexCom_Split_Pos.list[i]))
        attach.df$Genes <- as.factor(attach.df$Genes)
        CCMar_SexCom_Split_Pos_CT.df <- full_join(CCMar_SexCom_Split_Pos_CT.df, 
                                                  attach.df,
                                                  by="Genes")
       
      }
    }
    df_Count <- CCMar_SexCom_Split_Pos_CT.df
    
  }else{# mode2 Count genes by combined all sub type
      df <- do.call(rbind, lapply(CCMar_SexCom_Split_Pos.list, as.data.frame)) 
      colnames(df )[1]="Genes"
      df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
  }
  return(df_Count)
}



## Function of count gene number ##
  CountGene<- function(CCMar_SexCom_Split.list){
    df <- do.call(rbind, lapply(CCMar_SexCom_Split.list, as.data.frame))
    colnames(df )[1]="Genes"
    df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
    return(df_Count)
  }
# CCMar_SexCom_Count.df <- CountGene(CCMar_SexCom_Split.list)

# #Ori
# CCMar_SexCom_Split.df <- do.call(rbind, lapply(CCMar_SexCom_Split.list, as.data.frame)) 
# colnames(CCMar_SexCom_Split.df )[1]="Genes"
# CCMar_SexCom_Split_Count.df <- as.data.frame(table(CCMar_SexCom_Split.df$Genes)) %>% arrange(.,desc(Freq))
