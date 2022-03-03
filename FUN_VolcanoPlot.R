VolcanoPlot <- function(Tcell.Cachexia.Marker.F.S, Tcell.Cachexia.Marker.F.S_Pos_List, Tcell.Cachexia.Marker.F.S_Neg_List,
                        color = c(red = "#ef476f",gray = "gray",blue = "#0077b6"),
                        Xintercept = c(-1, 1),Yintercept = -log10(0.05),
                        ShowGeneNum = 5){
#
library(ggplot2)
library(cowplot)
##-------------- Volcano Plot --------------##
Tcell.Cachexia.Marker.F.S2 <- data.frame(row.names(Tcell.Cachexia.Marker.F.S),Tcell.Cachexia.Marker.F.S)
colnames(Tcell.Cachexia.Marker.F.S2)[[1]] <- c("Gene")

Tcell.Cachexia.Marker.F.S2$p_val <- Tcell.Cachexia.Marker.F.S2$p_val+1.0e-300

Tcell.Cachexia.Marker.F.S2$color <- ifelse(Tcell.Cachexia.Marker.F.S2$p_val<0.05 & abs(Tcell.Cachexia.Marker.F.S2$avg_log2FC)>= 1,ifelse(Tcell.Cachexia.Marker.F.S2$avg_log2FC > 1,'red','blue'),'gray')
# color <- c(red = "red",gray = "gray",blue = "blue")
# color <- c(red = "#ef476f",gray = "gray",blue = "#0077b6")

# redefine levels:
# Tcell.Cachexia.Marker.F.S2$genelabels <- factor(Tcell.Cachexia.Marker.F.S2$Gene, levels = c(Tcell.Cachexia.Marker.F.S_Pos_List,Tcell.Cachexia.Marker.F.S_Neg_List))

if (length(Tcell.Cachexia.Marker.F.S_Pos_List)>=ShowGeneNum && length(Tcell.Cachexia.Marker.F.S_Neg_List)>=ShowGeneNum) {
  Tcell.Cachexia.Marker.F.S2$genelabels <- factor(Tcell.Cachexia.Marker.F.S2$Gene, levels = c(Tcell.Cachexia.Marker.F.S_Pos_List[1:ShowGeneNum],Tcell.Cachexia.Marker.F.S_Neg_List[(length(Tcell.Cachexia.Marker.F.S_Neg_List)-(ShowGeneNum-1)):length(Tcell.Cachexia.Marker.F.S_Neg_List)]))
}  else {
  Tcell.Cachexia.Marker.F.S2$genelabels <- factor(Tcell.Cachexia.Marker.F.S2$Gene, levels = c(Tcell.Cachexia.Marker.F.S_Pos_List[1:length(Tcell.Cachexia.Marker.F.S_Pos_List)],Tcell.Cachexia.Marker.F.S_Neg_List[(length(Tcell.Cachexia.Marker.F.S_Neg_List)-(length(Tcell.Cachexia.Marker.F.S_Neg_List)-1)):length(Tcell.Cachexia.Marker.F.S_Neg_List)]))
  
}

library(ggrepel)
VolcanoPlot <- ggplot(Tcell.Cachexia.Marker.F.S2, aes(avg_log2FC, -log10(p_val), label = genelabels, col = color)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (p-value)") +
  geom_hline(yintercept = Yintercept, lty=8,col="black",lwd=0.8) +
  geom_vline(xintercept = Xintercept, lty=8,col="black",lwd=0.8) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        text = element_text(size = 15)) + 
  geom_point() +
  geom_text_repel(col = "#14213d", na.rm = TRUE,size = 6, box.padding = unit(0.45, "lines"), hjust = 1)+ 
  #geom_label(nudge_y = 2, alpha = 0.5)+
  theme(aspect.ratio=1)

VolcanoPlot_2 <- VolcanoPlot + theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
# https://www.cnblogs.com/liujiaxin2018/p/14257944.html #f!??????g40


# UMAP3 <- FeaturePlot(PBMC.combined, features = Tcell.Cachexia.Marker.F.S_Pos_List[1], split.by = SplitBy, max.cutoff = 3,
#                      cols = c("grey","#de3767", "red"), ncol = 2)
# UMAP4 <- FeaturePlot(PBMC.combined, features = Tcell.Cachexia.Marker.F.S_Neg_List[length(Tcell.Cachexia.Marker.F.S_Neg_List)], split.by = SplitBy, max.cutoff = 3,
#                      cols = c("grey", "blue"), ncol = 2)

return(VolcanoPlot_2)
}



