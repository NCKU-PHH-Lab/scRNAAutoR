# ## Reference: https://medium.com/datainpoint/r-essentials-writing-functions-f32d3c5cfc01
# # How to define a function
# FUNCTION_NAME <- function(INPUT1, INPUT2, ..., PARAM1, PARAM2, ...) {
#   # BODY
#   return(OUTPUT)
# }

##### Abbreviation Note #####
# Cachexia.Marker: CCMarker


##### Try #####
# ident1 <- "T_EO_Female"
# ident2 <-"T_LO_Female"


#####----------------------------- Function -----------------------------#####

Find_Markers <- function(PBMC.combined, ident1, ident2, CellType ,Path = getwd(), 
                         log2FC=1, Pval=0.05,ResultFolder = "FindMarkers",ProjectTitle="Tar"
                         ) {
  ## Find Cachexia marker in T Cell
  
  # Set empty data frame
  Tcell.CCMarker.F.S <- data.frame(p_val="", avg_log2FC="", pct.1="", pct.2="",
                                   p_val_adj="")
  Tcell.CCMarker.F.S <- Tcell.CCMarker.F.S[-1,]
  Tcell.CCMarker.F.All <-  Tcell.CCMarker.F.S
  Tcell.CCMarker.F.S_Order <- Tcell.CCMarker.F.S
  Tcell.CCMarker.F.S_Pos <- Tcell.CCMarker.F.S
  Tcell.CCMarker.F.S_Neg <- Tcell.CCMarker.F.S  
  Tcell.CCMarker.F.S_Pos_List <- ""
  Tcell.CCMarker.F.S_Neg_List <- ""
  
  try({
  set.seed(1) # Fix the seed
  Tcell.CCMarker.F.S <- FindMarkers(PBMC.combined, 
                                           ident.1 = ident1, ident.2 = ident2, 
                                           verbose = FALSE)
  
  
  Tcell.CCMarker.F.All <- FindMarkers(PBMC.combined,
                                             ident.1 = ident1, ident.2 = ident2,
                                             logfc.threshold = 0,min.pct=0,
                                             verbose = FALSE)
  
  Tcell.CCMarker.F.All <- data.frame(Genes = row.names(Tcell.CCMarker.F.All),Tcell.CCMarker.F.All)
  })
  
  write.table(Tcell.CCMarker.F.All, file=paste0(Save.Path,"/", ResultFolder,"/", CellType,"_FindMarkers.tsv"),
              sep="\t", row.names= F, quote = FALSE)
  try({
  Tcell.CCMarker.F.S_Order <- Tcell.CCMarker.F.S[order(Tcell.CCMarker.F.S$avg_log2FC,decreasing=TRUE, na.last=FALSE),]
  
  Tcell.CCMarker.F.S_Pos <- Tcell.CCMarker.F.S_Order[Tcell.CCMarker.F.S_Order$avg_log2FC>=log2FC 
                                                                   & Tcell.CCMarker.F.S_Order$p_val < Pval,]
  Tcell.CCMarker.F.S_Pos_List <- row.names(Tcell.CCMarker.F.S_Pos)
  
  Tcell.CCMarker.F.S_Neg <- Tcell.CCMarker.F.S_Order[Tcell.CCMarker.F.S_Order$avg_log2FC<=-log2FC 
                                                                   & Tcell.CCMarker.F.S_Order$p_val < Pval,]
  Tcell.CCMarker.F.S_Neg_List <- row.names(Tcell.CCMarker.F.S_Neg)
  
  # #
  # library(ggplot2)
  # library(cowplot)
  # 
  # 
  # # UMAP1 <-FeaturePlot(PBMC.combined, features = Tcell.CCMarker.F.S_Pos_List[1:2], split.by = SplitBy, max.cutoff = 3,
  # #                     cols = c("grey","#de3767", "red"))
  # # 
  # # UMAP2 <-FeaturePlot(PBMC.combined, features = Tcell.CCMarker.F.S_Neg_List[(length(Tcell.CCMarker.F.S_Neg_List)-1):length(Tcell.CCMarker.F.S_Neg_List)], split.by = SplitBy, max.cutoff = 3,
  # #                     cols = c("grey", "blue"))
  # # #
  # # VlnPlot1 <- VlnPlot(PBMC.combined, features = Tcell.CCMarker.F.S_Pos_List[1:2], split.by = SplitBy, group.by = "celltype",
  # #                     pt.size = 0, combine = FALSE, cols=c( "#d60b66","#5129f0","#f06413","#0f9ad6"))
  # # VlnPlot1_2 <- wrap_plots(plots = VlnPlot1, ncol = 2)
  # # 
  # # 
  # # VlnPlot2 <- VlnPlot(PBMC.combined, 
  # #                     features = Tcell.CCMarker.F.S_Neg_List[(length(Tcell.CCMarker.F.S_Neg_List)-1):length(Tcell.CCMarker.F.S_Neg_List)], 
  # #                     split.by = SplitBy, group.by = "celltype",
  # #                     pt.size = 0, combine = FALSE, cols=c( "#d60b66","#5129f0","#f06413","#0f9ad6"))
  # # VlnPlot2_2 <- wrap_plots(plots = VlnPlot2, ncol = 2)
  # # 
  # ##-------------- Volcano Plot --------------##
  # Tcell.CCMarker.F.S2 <- data.frame(row.names(Tcell.CCMarker.F.S),Tcell.CCMarker.F.S)
  # colnames(Tcell.CCMarker.F.S2)[[1]] <- c("Gene")
  # 
  # Tcell.CCMarker.F.S2$p_val <- Tcell.CCMarker.F.S2$p_val+1.0e-300
  # 
  # Tcell.CCMarker.F.S2$color <- ifelse(Tcell.CCMarker.F.S2$p_val<0.05 & abs(Tcell.CCMarker.F.S2$avg_log2FC)>= 1,ifelse(Tcell.CCMarker.F.S2$avg_log2FC > 1,'red','blue'),'gray')
  # # color <- c(red = "red",gray = "gray",blue = "blue")
  # # color <- c(red = "#ef476f",gray = "gray",blue = "#0077b6")
  # 
  # # redefine levels:
  # # Tcell.CCMarker.F.S2$genelabels <- factor(Tcell.CCMarker.F.S2$Gene, levels = c(Tcell.CCMarker.F.S_Pos_List,Tcell.CCMarker.F.S_Neg_List))
  # 
  # if (length(Tcell.CCMarker.F.S_Pos_List)>=5 && length(Tcell.CCMarker.F.S_Neg_List)>=5) {
  #   Tcell.CCMarker.F.S2$genelabels <- factor(Tcell.CCMarker.F.S2$Gene, levels = c(Tcell.CCMarker.F.S_Pos_List[1:5],Tcell.CCMarker.F.S_Neg_List[(length(Tcell.CCMarker.F.S_Neg_List)-4):length(Tcell.CCMarker.F.S_Neg_List)]))
  # }  else if (length(Tcell.CCMarker.F.S_Pos_List)>=2 && length(Tcell.CCMarker.F.S_Neg_List)>=2) {
  #   Tcell.CCMarker.F.S2$genelabels <- factor(Tcell.CCMarker.F.S2$Gene, levels = c(Tcell.CCMarker.F.S_Pos_List[1:2],Tcell.CCMarker.F.S_Neg_List[(length(Tcell.CCMarker.F.S_Neg_List)-1):length(Tcell.CCMarker.F.S_Neg_List)]))
  #   
  # }
  # else {
  #   Tcell.CCMarker.F.S2$genelabels <- factor(Tcell.CCMarker.F.S2$Gene, levels = c(Tcell.CCMarker.F.S_Pos_List[1],Tcell.CCMarker.F.S_Neg_List[length(Tcell.CCMarker.F.S_Neg_List)]))
  # }
  # 
  # library(ggrepel)
  # VolcanoPlot <- ggplot(Tcell.CCMarker.F.S2, aes(avg_log2FC, -log10(p_val), label = genelabels, col = color)) +
  #   geom_point(size = 3) +
  #   theme_bw() +
  #   scale_color_manual(values = color) +
  #   labs(x="log2 (fold change)",y="-log10 (p-value)") +
  #   geom_hline(yintercept = -log10(0.05), lty=8,col="black",lwd=0.8) +
  #   geom_vline(xintercept = c(-1, 1), lty=8,col="black",lwd=0.8) +
  #   theme(legend.position = "none",
  #         panel.grid=element_blank(),
  #         axis.title = element_text(size = 16),
  #         axis.text = element_text(size = 14),
  #         text = element_text(size = 15)) + 
  #   geom_point() +
  #   geom_text_repel(col = "#14213d", na.rm = TRUE,size = 6, box.padding = unit(0.45, "lines"), hjust = 1)+ 
  #   #geom_label(nudge_y = 2, alpha = 0.5)+
  #   theme(aspect.ratio=1)
  # 
  # VolcanoPlot_2 <- VolcanoPlot + theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  # # https://www.cnblogs.com/liujiaxin2018/p/14257944.html #f!??????g40
  # 
  # 
  # # UMAP3 <- FeaturePlot(PBMC.combined, features = Tcell.CCMarker.F.S_Pos_List[1], split.by = SplitBy, max.cutoff = 3,
  # #                      cols = c("grey","#de3767", "red"), ncol = 2)
  # # UMAP4 <- FeaturePlot(PBMC.combined, features = Tcell.CCMarker.F.S_Neg_List[length(Tcell.CCMarker.F.S_Neg_List)], split.by = SplitBy, max.cutoff = 3,
  # #                      cols = c("grey", "blue"), ncol = 2)
  # 
  # 
  # 
  # 
  # ##------------------------ Output ------------------------ ##
  # # OUTPUT <- list(Tcell.CCMarker.F.S, Tcell.CCMarker.F.All, 
  # #                Tcell.CCMarker.F.S_Pos, Tcell.CCMarker.F.S_Pos_List,
  # #                Tcell.CCMarker.F.S_Neg, Tcell.CCMarker.F.S_Neg_List,
  # #                UMAP1, UMAP2, VlnPlot1_2, VlnPlot2_2, VolcanoPlot_2)
  # # names(OUTPUT) <- c('CCMarker.S','CCMarker.All',
  # #                    'CCMarker.S_Pos', 'CCMarker.S_Pos_List',
  # #                    'CCMarker.S_Neg', 'CCMarker.S_Neg_List',
  # #                    "UMAP1","UMAP2","VlnPlot1","VlnPlot2","VolcanoPlot")
  })
  
  OUTPUT <- list(Tcell.CCMarker.F.S, Tcell.CCMarker.F.All, 
                 Tcell.CCMarker.F.S_Pos, Tcell.CCMarker.F.S_Pos_List,
                 Tcell.CCMarker.F.S_Neg, Tcell.CCMarker.F.S_Neg_List)
  names(OUTPUT) <- c(paste0(ProjectTitle,'Marker.S'), paste0(ProjectTitle,'Marker.All'),
                     paste0(ProjectTitle,'Marker.S_Pos'), paste0(ProjectTitle,'Marker.S_Pos_List'),
                     paste0(ProjectTitle,'Marker.S_Neg'), paste0(ProjectTitle,'Marker.S_Neg_List'))
  
  
  return(OUTPUT)
}

