# ##### RNA-seq analysis in R #####
# # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
# library(fgsea)
# gseaDat <- Cachexia.Marker.Male[["Cachexia.Marker.Male.Mac"]][["Cachexia.Marker.All"]]
# gseaDat <- data.frame(row.names(gseaDat),gseaDat)
# colnames(gseaDat)[[1]] <- c("Gene")
# ranks <- gseaDat$avg_log2FC
# names(ranks) <- gseaDat$Gene
# head(ranks)
# 
# barplot(sort(ranks, decreasing = T))
# 
# # Geneset from GSEA
# Pathway.all <- read.delim(paste0(PathName,"/Pathway.all.v7.4.symbols.gmt"),header = F)

FUN_GSEA_LargeGeneSet <- function(ranks,Pathway.all,TopNum) {
  library(fgsea)  

    Pathway.all.list <- list()
        for (i in c(1:length(Pathway.all[,1]))) {
          try({
          
          Pathway.all.list.ori <- as.data.frame(t(Pathway.all[i,3:length(Pathway.all[i,])]))
          colnames(Pathway.all.list.ori)[[1]] <- c("Gene")
          
          ## tryCatch(
          ##   {
          Pathway.all.list.ori <- na.omit(Pathway.all.list.ori)
          ###Pathway.all.list.ori <- HSsymbol2MMsymbol(Pathway.all.list.ori,"Gene")
         
          # Delete NA(or 0)
          Pathway.all.list.ori <- Pathway.all.list.ori[Pathway.all.list.ori!=0,]
          #Bug# Pathway.all.list.ori <- Pathway.all.list.ori[-which(Pathway.all.list.ori$MM.symbol==0),]
          #Bug# Pathway.all.list.ori <- Pathway.all.list.ori[!which(Pathway.all.list.ori$MM.symbol==0),]
          #NoUse# Pathway.all.list[[i]] <- na.omit(Pathway.all.list[[i]])
          #Error# Pathway.all.list[[i]] <- Pathway.all.list[!is.na(Pathway.all.list)]
          
          Pathway.all.list.ori <- unique(Pathway.all.list.ori)
          ###Pathway.all.list.ori <- unique(Pathway.all.list.ori$MM.symbol)
          Pathway.all.list[[i]] <- as.character(Pathway.all.list.ori)
          # Pathway.all.list[[i]] <- as.character(Pathway.all[i,3:length(Pathway.all[i,])])  
          
          
          names(Pathway.all.list)[[i]] <- Pathway.all[i,1]
          rm(Pathway.all.list.ori)
          ##  },
          ## error=function(e) {
          ## Pathway.all.list <- Pathway.all.list
          ## })
          })
        }
    
    # load(paste0(PathName,"/Full_annotation.RData"))
    pathwaysH <- Pathway.all.list 
    
    fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500)
    # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    
    # TOP N
    library(magrittr)
    library(dplyr)
    topUp <- fgseaRes %>% 
      dplyr::filter(ES > 0) %>% 
      top_n(TopNum, wt=-padj)
    topDown <- fgseaRes %>% 
      dplyr::filter(ES < 0) %>% 
      top_n(TopNum, wt=-padj)
    topPathways <- bind_rows(topUp, topDown) %>% 
      arrange(-ES)
  #####------------------------ Output ------------------------ #####
  OUTPUT <- list(Pathway.all.list,pathwaysH,fgseaRes,topPathways)
  names(OUTPUT) <- c("Pathway.all.list","pathwaysH","fgseaRes","topPathways")
  
  rm(pathwaysH,fgseaRes,Pathway.all.list)
return(OUTPUT)
}


# GSEAOutput <- FUN_GSEA_LargeGeneSet(gseaDat, ranks,Pathway.all)
# 
# fgseaRes <- GSEAOutput[["fgseaRes"]]
# head(fgseaRes[order(padj, -abs(NES)), ], n=10)
# plot.new()
# plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
# 
# library(magrittr)
# library(dplyr)
# topUp <- fgseaRes %>% 
#   filter(ES > 0) %>% 
#   top_n(10, wt=-padj)
# topDown <- fgseaRes %>% 
#   filter(ES < 0) %>% 
#   top_n(10, wt=-padj)
# topPathways <- bind_rows(topUp, topDown) %>% 
#   arrange(-ES)
# plot.new()
# plotGseaTable(pathwaysH[topPathways$pathway], 
#               ranks, 
#               fgseaRes, 
#               gseaParam = 0.5)


