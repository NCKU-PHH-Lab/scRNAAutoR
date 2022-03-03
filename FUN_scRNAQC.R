scRNAQC <- function(PBMC.combined, nUMIFilter= 500, nGeneFilter = 250 ,
                    logGPUFilter= 0.8 , mitRFilter=0.1 ,
                    PtSize=0,
                    AddMitInf = "Yes", # Add mitochondria information
                    CheckOnly="No", # CheckOnly = "Yes": Just plot the Fig
                    FileName = "QC"){   

    ## https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
    ## QC for scRNA Data
    
    ##### QC in Seruat #####
    PBMC.combined_QC <- PBMC.combined
    # QC and selecting cells for further analysis
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    ## PBMC.combined_QC[["percent.mt"]] <- PercentageFeatureSet(PBMC.combined_QC, pattern = "^MT-") # For Human
    ## PBMC.combined_QC[["percent.mt"]] <- PercentageFeatureSet(PBMC.combined_QC, pattern = "^mt-") # For Mouse
    if(AddMitInf == "Yes"){
        PBMC.combined_QC <- scRNAMit(PBMC.combined_QC)
      }else{
        PBMC.combined_QC <- PBMC.combined_QC
      }
        # Visualize QC metrics as a violin plot
        # VlnPlot(PBMC.combined_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # from Seruat
        VlnPlot.PBMC <- VlnPlot(PBMC.combined_QC, features = c("nGene", "nUMI", "mitoRatio"), ncol = 2,pt.size= PtSize)
        FSplot1.PBMC <- FeatureScatter(PBMC.combined_QC, feature1 = "nUMI", feature2 = "nGene")
        VlnPlot.PBMC + FSplot1.PBMC
        
        
        #########################################################################################################
        ##### QC in Harvard Chan Bioinformatics Core (HBC) #####
        ## https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
        # Add number of genes per UMI for each cell to metadata
        PBMC.combined_QC$log10GenesPerUMI <- log10(PBMC.combined_QC$nGene) / log10(PBMC.combined_QC$nUMI)
        
        metadata <- PBMC.combined_QC@meta.data
        
        ## Cell counts
        # Visualize the number of cell counts per sample
        metadata %>% 
          ggplot(aes(x=sample, fill=sample)) + 
          geom_bar() +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),aspect.ratio=1) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells") -> NCells.BarPlot
        
        NCells.BarPlot
        
        
        ## UMI counts (transcripts) per cell
        # Visualize the number UMIs/transcripts per cell
        metadata %>% 
          ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
          geom_density(alpha = 0.2) + 
          scale_x_log10() + 
          theme_classic() +
          ylab("Cell density") +
          geom_vline(xintercept = nUMIFilter) -> Cell_density.UMIs.Plot
        
        Cell_density.UMIs.Plot
        
        
        ## Genes detected per cell
        # Visualize the distribution of genes detected per cell via histogram
        metadata %>% 
          ggplot(aes(color=sample, x=nGene, fill= sample)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() +
          ylab("Cell density") +
          geom_vline(xintercept = nGeneFilter) ->  Cell_density.NGenes.Plot
        
        Cell_density.NGenes.Plot
        
        # Visualize the distribution of genes detected per cell via boxplot
        metadata %>% 
          ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
          geom_boxplot() + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),aspect.ratio=1) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes") ->  NCvsNG.BarPlot
        
        NCvsNG.BarPlot
        
        
        ## UMIs vs. genes detected
        
        # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
        metadata %>% 
          ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black") +
          stat_smooth(method=lm) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          geom_vline(xintercept = nUMIFilter) +
          geom_hline(yintercept = nGeneFilter) +
          facet_wrap(~sample) -> FSplot
        
        FSplot
        
        
        ## Mitochondrial counts ratio
        
        # Visualize the distribution of mitochondrial gene expression detected per cell
        metadata %>% 
          ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
          geom_density(alpha = 0.2) + 
          #scale_x_log10() + 
          theme_classic() +
          ylab("Cell density") +
          geom_vline(xintercept = 0.1) ->  Cell_density.mitR.Plot
        
        Cell_density.mitR.Plot
        
        
        ## Complexity
        # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
        metadata %>%
          ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
          geom_density(alpha = 0.2) +
          theme_classic() +
          geom_vline(xintercept = logGPUFilter) -> Cell_density.log10GenesPerUMI.Plot
        
        Cell_density.log10GenesPerUMI.Plot
        
        # Cell_density
        Cell_density.UMIs.Plot + Cell_density.NGenes.Plot + 
          Cell_density.mitR.Plot + Cell_density.log10GenesPerUMI.Plot
        # BarPlot
        NCells.BarPlot + NCvsNG.BarPlot
        # FeatureScatter
        FSplot 
    
        # https://stackoverflow.com/questions/19288101/r-pdf-usage-inside-a-function/19288874
        pdf(
          file = paste0(getwd(),"/",FileName,".pdf"),
          width = 10,  height = 7
        )
        
        print(VlnPlot.PBMC + FSplot1.PBMC)
        # Cell_density
        print(Cell_density.UMIs.Plot + Cell_density.NGenes.Plot + Cell_density.mitR.Plot + Cell_density.log10GenesPerUMI.Plot)
        # BarPlot
        print(NCells.BarPlot + NCvsNG.BarPlot)
        # FeatureScatter
        print(FSplot)
        
        dev.off()
        # graphics.off()
        
      if(CheckOnly == "No"){
        ##### Filtering #####
        ## Cell-level filtering
        # Filter out low quality reads using selected thresholds - these will change with experiment
        PBMC.combined_QC_Filter <- subset(x = PBMC.combined_QC, 
                                          subset= (nUMI >= nUMIFilter) & 
                                            (nGene >= nGeneFilter) & 
                                            (log10GenesPerUMI > logGPUFilter) & 
                                            (mitoRatio < mitRFilter))
        
        # ## Gene-level filtering
        # # Output a logical vector for every gene on whether the more than zero counts per cell
        # # Extract counts
        # # counts <- GetAssayData(object = PBMC.combined_QC_Filter, slot = "counts")
        # counts <- PBMC.combined_QC_Filter@assays[["RNA"]]@counts
        # 
        # 
        # # Output a logical vector for every gene on whether the more than zero counts per cell
        # nonzero <- counts > 0
        # 
        # # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
        # keep_genes <- Matrix::rowSums(nonzero) >= 10
        # 
        # # Only keeping those genes expressed in more than 10 cells
        # filtered_counts <- counts[keep_genes, ]
        # 
        # # Reassign to filtered Seurat object
        # PBMC.combined_QC_Filter2 <- CreateSeuratObject(filtered_counts, meta.data = PBMC.combined_QC_Filter@meta.data)
        # scRNAQC(PBMC.combined_QC_Filter2,AddMitInf = "No",FileName = paste0(FileName,"_CheckTry"), CheckOnly="Yes")
        scRNAQC(PBMC.combined_QC_Filter,AddMitInf = "No",FileName = paste0(FileName,"_Check"), CheckOnly="Yes")
      }else{
        PBMC.combined_QC_Filter <- PBMC.combined_QC
      }


return(PBMC.combined_QC_Filter)
}



##### QC Analysis for single Sample #####
# # TN138
# PBMC.TN138.combined_QC <- PBMC.TN138
# # QC and selecting cells for further analysis
# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# PBMC.TN138.combined_QC[["percent.mt"]] <- PercentageFeatureSet(PBMC.TN138.combined_QC, pattern = "^MT-")
# PBMC.TN138.combined_QC[["percent.mt"]] <- PercentageFeatureSet(PBMC.TN138.combined_QC, pattern = "^mt-") # OK, and same as mitoRatio
# 
# # Visualize QC metrics as a violin plot
# v1 <- VlnPlot(PBMC.TN138.combined_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# p1 <- FeatureScatter(PBMC.TN138, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p2 <- FeatureScatter(PBMC.TN138.combined_QC, feature1 = "nCount_RNA", feature2 = "percent.mt")
# p3 <- p1+p2
# 
# 
# # Test by funciton
# PBMC.TN138.combined_QC2 <- scRNAMit(PBMC.TN138)
# v1T <- VlnPlot(PBMC.TN138.combined_QC2, features = c("nGene", "nUMI", "mitoRatio"), ncol = 3)
# v1/v1T
# p1T <- FeatureScatter(PBMC.TN138.combined_QC2, feature1 = "nUMI", feature2 = "nGene")
# p2T <- FeatureScatter(PBMC.TN138.combined_QC2, feature1 = "nUMI", feature2 = "mitoRatio")
# p3T <- p1T+p2T
# p3/p3T

