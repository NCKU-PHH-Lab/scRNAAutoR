GSEA_ggplot = function( GSEA_Large.Sum.TOP,
                        GSEA_Color = list(high = "#ef476f",mid = "white",low = "#0077b6"),
                        NES_Th = 1, # Th:Threshold
                        padj_Th = 0.05
){

  ##### Bubble plot #####
  GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > NES_Th,]
  GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < padj_Th,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]
  library(ggplot2)
  library(scales)
  
  BBPlot_Ori <- ggplot(GSEA_Large.Sum.TOP.S,aes(x = PhenoType, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() + scale_colour_gradient2(low = GSEA_Color[["low"]], mid = GSEA_Color[["mid"]], high = GSEA_Color[["high"]],  # low = "#04873f", mid = "white", high = "#e3672d"
                                          guide = "colourbar",midpoint = 0)+
    theme(axis.text.x = element_text(face="bold",  size = 12 ,angle = 90, hjust = 1, vjust = .5))
  # BBPlot_Ori <- ggplot(GSEA_Large.Sum.TOP.S,aes(x = GSEA_Large.Sum.TOP.S[,"PhenoType"], y = GSEA_Large.Sum.TOP.S[,"pathway"], 
  #                                               color = GSEA_Large.Sum.TOP.S[,"NES"], size = -log10(GSEA_Large.Sum.TOP.S[,"padj"]))) + 
  #   geom_point() + scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f",  # low = "#04873f", mid = "white", high = "#e3672d"
  #                                         guide = "colourbar",midpoint = 0)+            
  #   theme(axis.text.x = element_text(face="bold",  size = 12 ,angle = 90, hjust = 1, vjust = .5))
  
  BBPlot_Ori
  BBPlot <- BBPlot_Ori + geom_point() + 
    theme(legend.position = "bottom") + theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BBPlot
  BBPlot <- BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                      XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
  BBPlot
  
  BBPlot2 <- BBPlot +theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank())
  BBPlot2
  # print(BBPlot)
  # print(BBPlot2)
  
  
  
  ##### Clustering #####
  #GSEA.df <- GSEA_Large.df
  GSEA.df <- GSEA_Large.Sum.TOP.S
  
  
  # GSEA_Color.lt = list(red = "#ef476f",gray = "gray",blue = "#0077b6")
  GSEA_Color.lt = list(red = "#04873f",gray = "gray",blue = "#e3672d")
  
  GSEA.df$PhenoType <- factor(GSEA.df$PhenoType,
                              levels = sort(unique(as.character(GSEA.df$PhenoType))))
  
  # BBPlot <- ggplot(GSEA.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + geom_point() + 
  #   scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
  #                          guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ scale_x_discrete()+ theme_bw()+
  #   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
  # 
  # 
  # 
  # BBPlot
  # 
  # BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
  #                           XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
  
  
  df1 <- reshape2::dcast(GSEA.df,PhenoType~pathway,value.var = "NES")
  rownames(df1)<-df1$id
  
  df1.1<-df1[,2:ncol(df1)]
  df1.1[is.na(df1.1)] <- 0
  
  
  df1.1.clust.Pheno<-hclust(dist(df1.1))
  df1.1.clust.Pathway<-hclust(dist(t(df1.1)))
  
  
  library(ggtree)
  library(dplyr)
  PhenoType_Order <- data.frame(No=row.names(df1),PhenoType=df1[,1]) 
  PhenoType_Order$No <- as.numeric(PhenoType_Order$No)
  PhenoType_Treeclust_Order <- data.frame(No=df1.1.clust.Pheno[["order"]])
  PhenoType_Treeclust_Order <- left_join(PhenoType_Treeclust_Order,PhenoType_Order) 
  df1.1.clust.Pheno[["order"]] <- PhenoType_Treeclust_Order$PhenoType
  df1.1.clust.Pheno[["labels"]] <- PhenoType_Treeclust_Order$PhenoType
  
  p2 <- ggtree(df1.1.clust.Pheno)
  p2+
    geom_tiplab()+
    xlim(NA,7)
  p2.2<-p2+
    geom_tiplab()+
    xlim(NA,7)+
    #geom_tiplab(angle=90)+
    #theme_tree2()+
    layout_dendrogram()
  p2.2
  
  p3<-ggtree(df1.1.clust.Pathway)
  p3+
    geom_tiplab()+
    xlim(NA,7)
  ## Plot
  library(aplot)
  BBPlot_Cluster<- BBPlot %>%
    insert_left(p3,width = 0.2)
  
  BBPlot_Cluster
  
  BBPlotB <- BBPlot %>% 
    BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                   XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=2.5, XaThick=0.6, YaThick=0.6,
                   LegTextSize = 12,LegTitleSize=13)
  BBPlotB1 <- BBPlotB %>%
    insert_left(p3,width = 0.5)
  BBPlotB1
  
  ## Dot size
  # https://community.rstudio.com/t/scale-geom-point-size-to-increase-size-based-on-distance-from-zero/20244/2
  
  
  
  # BBPlotB2 <-BBPlotB %>%
  #   insert_left(p3,width = 0.2)%>%
  #   insert_top(p2+layout_dendrogram(),height = 0.2)
  # BBPlotB2

  GSEA_BBPlot.lt <- list(GSEA_Large.Sum.TOP.S, BBPlot_Ori, BBPlot, BBPlot2, BBPlotB1,p3)
  names(GSEA_BBPlot.lt) <- c("GSEA_TOP.df", "BBPlot_Ori", "BBPlot", "BBPlot2", "BBPlotB1", "Y_Order")
  # pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_All.pdf"),width = 17, height = 10 )
  #   BBPlot_Ori
  #   BBPlot
  #   BBPlot2
  #   BBPlotB1
  #   #BBPlotB2
  # dev.off()
  
  return(GSEA_BBPlot.lt)
}