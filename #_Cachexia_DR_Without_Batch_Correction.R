## 00_Cachexia_Dimension_Reduction_Without_Batch_Correction.R

##### Version information ######
  # _                           
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          1.1                         
  # year           2021                        
  # month          08                          
  # day            10                          
  # svn rev        80725                       
  # language       R                           
  # version.string R version 4.1.1 (2021-08-10)
  # nickname       Kick Things  

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting  ##### 
  Version = "20220225_BatchEffect"
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)
  
##### Load Packages  ##### 
  library(Seurat)
  library(ggplot2)
  library("dplyr")

##### Function setting  ##### 
  ## Call function
  source("FUN_Beautify_ggplot.R")

#####-------------- SC --------------#####
  ##### Load datasets  #####
    SC.data.TN136 <- Read10X(data.dir = paste0(getwd(),"/TN136/monocle/outs/filtered_gene_bc_matrices/mm10"))
    SC.TN136 <- CreateSeuratObject(counts = SC.data.TN136, project = "EO.M", min.cells = 3, min.features = 200)
    SC.TN136
    SC.TN136@meta.data[["sample"]] <- rep(c("EO.M"), times=length(SC.TN136@meta.data[["orig.ident"]]))
    SC.TN136@meta.data[["ID"]] <- rep(c("TN136"), times=length(SC.TN136@meta.data[["orig.ident"]]))
    SC.TN136@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN136@meta.data[["orig.ident"]]))  #EO: Early_Onset
    SC.TN136@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN136@meta.data[["orig.ident"]]))
    
    SC.data.TN137 <- Read10X(data.dir = paste0(getwd(),"/TN137/monocle/outs/filtered_gene_bc_matrices/mm10"))
    SC.TN137 <- CreateSeuratObject(counts = SC.data.TN137, project = "LO.M", min.cells = 3, min.features = 200)
    SC.TN137
    SC.TN137@meta.data[["sample"]] <- rep(c("LO.M"), times=length(SC.TN137@meta.data[["orig.ident"]]))
    SC.TN137@meta.data[["ID"]] <- rep(c("TN137"), times=length(SC.TN137@meta.data[["orig.ident"]]))
    SC.TN137@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN137@meta.data[["orig.ident"]]))  #LO: Late_Onset
    SC.TN137@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN137@meta.data[["orig.ident"]]))
    
    SC.data.TN145 <- Read10X(data.dir = paste0(getwd(),"/TN145/monocle/outs/filtered_gene_bc_matrices/mm10"))
    SC.TN145 <- CreateSeuratObject(counts = SC.data.TN145, project = "LO.F", min.cells = 3, min.features = 200)
    SC.TN145
    SC.TN145@meta.data[["sample"]] <- rep(c("LO.F"), times=length(SC.TN145@meta.data[["orig.ident"]]))
    SC.TN145@meta.data[["ID"]] <- rep(c("TN145"), times=length(SC.TN145@meta.data[["orig.ident"]]))
    SC.TN145@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN145@meta.data[["orig.ident"]]))
    SC.TN145@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN145@meta.data[["orig.ident"]]))
    
    SC.data.TN147 <- Read10X(data.dir = paste0(getwd(),"/TN147/monocle/outs/filtered_gene_bc_matrices/mm10"))
    SC.TN147 <- CreateSeuratObject(counts = SC.data.TN147, project = "EO.F", min.cells = 3, min.features = 200)
    SC.TN147
    SC.TN147@meta.data[["sample"]] <- rep(c("EO.F"), times=length(SC.TN147@meta.data[["orig.ident"]]))
    SC.TN147@meta.data[["ID"]] <- rep(c("TN147"), times=length(SC.TN147@meta.data[["orig.ident"]]))
    SC.TN147@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN147@meta.data[["orig.ident"]]))
    SC.TN147@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN147@meta.data[["orig.ident"]]))

  ##### Merging More Than Two Seurat Objects #####  
    ## Merging More Than Two Seurat Objects
    ## https://satijalab.org/seurat/articles/merge_vignette.html  
    SC.combined <- merge(SC.TN136, y = c(SC.TN137, SC.TN145, SC.TN147), 
                         add.cell.ids = c("TN136", "TN137", "TN145","TN147"), 
                         project = "SC")
    SC.combined

  # ##### Perform integration #####  
  #   # this command creates an 'integrated' data assay
  #   SC.combined <- IntegrateData(anchorset = SC.combined)
  #   # specify that we will perform downstream analysis on the corrected data note that the
  #   # original unmodified data still resides in the 'RNA' assay
  #   DefaultAssay(SC.combined) <- "integrated"

  ## Perform an integrated analysis  
  # https://satijalab.org/seurat/articles/integration_introduction.html
    ## 
    # Run the standard workflow for visualization and clustering
    set.seed(1) # Fix the seed
    SC.combined <- ScaleData(SC.combined, verbose = FALSE)
    
    # RunPCA error
    # https://github.com/satijalab/seurat/issues/1788
    set.seed(1) # Fix the seed
    SC.combined <- FindVariableFeatures(object = SC.combined)
    set.seed(1) # Fix the seed
    SC.combined <- RunPCA(SC.combined, npcs = 30, verbose = FALSE)
    
    set.seed(1) # Fix the seed
    SC.combined <- RunUMAP(SC.combined, reduction = "pca", dims = 1:30)
    set.seed(1) # Fix the seed
    SC.combined <- FindNeighbors(SC.combined, reduction = "pca", dims = 1:30)
    set.seed(1) # Fix the seed
    SC.combined <- FindClusters(SC.combined, resolution = 0.5)
    
    # Visualization
    p1 <- DimPlot(SC.combined, reduction = "umap", group.by = "sample")
    p2 <- DimPlot(SC.combined, reduction = "umap", label = TRUE, repel = TRUE)
    p1 + p2  
    DimPlot(SC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
    

#####-------------- PBMC --------------#####
    ##### Load datasets  #####
      PBMC.data.TN138 <- Read10X(data.dir = paste0(getwd(),"/TN138/monocle/outs/filtered_gene_bc_matrices/mm10"))
      PBMC.TN138 <- CreateSeuratObject(counts = PBMC.data.TN138, project = "EO.M", min.cells = 3, min.features = 200)
      PBMC.TN138
      PBMC.TN138@meta.data[["sample"]] <- rep(c("EO.M"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))
      PBMC.TN138@meta.data[["ID"]] <- rep(c("TN138"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))
      PBMC.TN138@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))  #EO: Early_Onset
      PBMC.TN138@meta.data[["Sex"]] <- rep(c("Male"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))
      
      PBMC.data.TN139 <- Read10X(data.dir = paste0(getwd(),"/TN139/monocle/outs/filtered_gene_bc_matrices/mm10"))
      PBMC.TN139 <- CreateSeuratObject(counts = PBMC.data.TN139, project = "LO.M", min.cells = 3, min.features = 200)
      PBMC.TN139
      PBMC.TN139@meta.data[["sample"]] <- rep(c("LO.M"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))
      PBMC.TN139@meta.data[["ID"]] <- rep(c("TN139"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))
      PBMC.TN139@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))  #LO: Late_Onset
      PBMC.TN139@meta.data[["Sex"]] <- rep(c("Male"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))
      
      PBMC.data.TN146 <- Read10X(data.dir = paste0(getwd(),"/TN146/monocle/outs/filtered_gene_bc_matrices/mm10"))
      PBMC.TN146 <- CreateSeuratObject(counts = PBMC.data.TN146, project = "LO.F", min.cells = 3, min.features = 200)
      PBMC.TN146
      PBMC.TN146@meta.data[["sample"]] <- rep(c("LO.F"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
      PBMC.TN146@meta.data[["ID"]] <- rep(c("TN146"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
      PBMC.TN146@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
      PBMC.TN146@meta.data[["Sex"]] <- rep(c("Female"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
      
      PBMC.data.TN148 <- Read10X(data.dir = paste0(getwd(),"/TN148/monocle/outs/filtered_gene_bc_matrices/mm10"))
      PBMC.TN148 <- CreateSeuratObject(counts = PBMC.data.TN148, project = "EO.F", min.cells = 3, min.features = 200)
      PBMC.TN148
      PBMC.TN148@meta.data[["sample"]] <- rep(c("EO.F"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
      PBMC.TN148@meta.data[["ID"]] <- rep(c("TN148"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
      PBMC.TN148@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
      PBMC.TN148@meta.data[["Sex"]] <- rep(c("Female"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
    
    ##### Merging More Than Two Seurat Objects #####  
    ## Merging More Than Two Seurat Objects
    ## https://satijalab.org/seurat/articles/merge_vignette.html  
    PBMC.combined <- merge(PBMC.TN138, y = c(PBMC.TN139, PBMC.TN146, PBMC.TN148), 
                         add.cell.ids = c("TN138", "TN139", "TN146","TN148"), 
                         project = "SC")
    PBMC.combined
    
    # ##### Perform integration #####  
    #   # this command creates an 'integrated' data assay
    #   PBMC.combined <- IntegrateData(anchorset = PBMC.combined)
    #   # specify that we will perform downstream analysis on the corrected data note that the
    #   # original unmodified data still resides in the 'RNA' assay
    #   DefaultAssay(PBMC.combined) <- "integrated"
    
    ## Perform an integrated analysis  
    # https://satijalab.org/seurat/articles/integration_introduction.html
    ## 
    # Run the standard workflow for visualization and clustering
    set.seed(1) # Fix the seed
    PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE)
    
    # RunPCA error
    # https://github.com/satijalab/seurat/issues/1788
    set.seed(1) # Fix the seed
    PBMC.combined <- FindVariableFeatures(object = PBMC.combined)
    set.seed(1) # Fix the seed
    PBMC.combined <- RunPCA(PBMC.combined, npcs = 30, verbose = FALSE)
    
    set.seed(1) # Fix the seed
    PBMC.combined <- RunUMAP(PBMC.combined, reduction = "pca", dims = 1:30)
    set.seed(1) # Fix the seed
    PBMC.combined <- FindNeighbors(PBMC.combined, reduction = "pca", dims = 1:30)
    set.seed(1) # Fix the seed
    PBMC.combined <- FindClusters(PBMC.combined, resolution = 0.5)
    
    # Visualization
    p1 <- DimPlot(PBMC.combined, reduction = "umap", group.by = "sample")
    p2 <- DimPlot(PBMC.combined, reduction = "umap", label = TRUE, repel = TRUE)
    p1 + p2  
    DimPlot(PBMC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)

##### Export PDF #####
    pdf(
      file = paste0(setwd(getwd()),"/",Version,"/nlDR_Without_Batch_Correction.pdf"),
      width = 10,  height = 8
    )
    
    DimPlot(PBMC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.80, 0.15),AxisTitleSize=1.1)
    DimPlot(SC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.80, 0.15),AxisTitleSize=1.1)
    DimPlot(PBMC.combined, reduction = "umap", group.by = "sample") %>% 
      BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
    DimPlot(SC.combined, reduction = "umap", group.by = "sample") %>% 
      BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
    
    dev.off()
    
#####    
save.image(paste0(Save.Path,"/00_Cachexia_Dimension_Reduction_Without_Batch_Correction.RData"))     
    