##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          1.2                         
  # year           2021                        
  # month          11                          
  # day            01                          
  # svn rev        81115                       
  # language       R                           
  # version.string R version 4.1.2 (2021-11-01)
  # nickname       Bird Hippie 
  #--------------------------------#
  # R version 4.1.2
  # Seurat 4.0.2

##### Load Packages  ##### 
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  library(ggplot2)
  library("stringr")
  library("magrittr")
  library("dplyr")

  # In fuction
  library("AnnotationHub")
  library(ensembldb)
  
##### Function setting  ##### 
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_scRNAQC.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_Venn.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Beautify_UMAP.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

##### Current path and new folder setting  ##### 
  ## GSE103322 HNSC
  # Version = paste0(Sys.Date(),"_","HNSC")
  # Save.Path = paste0(getwd(),"/",Version)
  # dir.create(Save.Path)
  # 
  # FileName = paste0("GSE103322_HNSCC_all_data/HNSCC_all_data.txt")

  Version = paste0(Sys.Date(),"_","CC_PBMC")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)
  
  InputFolder = "Input_files_10x" 
  InputAnno = "PBMC_Ano.csv"
  ProjectName = "CC"
  
##### Load datasets  #####
  # ## GSE103322 HNSC
  # GeneExp.df <- read.table(InputFolder, header=T, row.names = 1, sep="\t")
  # #raw_counts <- GeneExp.df[6:nrow(GeneExp.df),]
  # scRNA.SeuObj <- CreateSeuratObject(counts = GeneExp.df[6:nrow(GeneExp.df),],
  #                              min.cells = 3, min.genes = 200, 
  #                              project = "HNSC")
  # scRNA.SeuObj <- NormalizeData(scRNA.SeuObj)
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  # scRNA.SeuObj@meta.data[["sample"]] <- GeneExp.df[5,] %>% as.character()
  
  # list_files.set <- list.files(InputFolder,full.names = T)
  # Nfiles = length(list_files.set)

  #### Cachexia
  list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  Feature.set <- colnames(list_files.df)[-1]
  
  scRNA_SeuObj.list <- list()
  for(i in 1:nrow(list_files.df)){
    Folder <- list_files.df$Folder[i]
    Data.dgCMatrix <- Read10X(data.dir = paste0(InputFolder,"/", Folder, "/monocle/outs/filtered_gene_bc_matrices/mm10"))
    Data.SeuObj <- CreateSeuratObject(counts = Data.dgCMatrix, project = ProjectName, min.cells = 3, min.features = 200)
    
    for (j in 1:(ncol(list_files.df)-1)) {
      Data.SeuObj@meta.data[[colnames(list_files.df)[j+1]]] <- rep(list_files.df[i,j+1], 
                                                                   times=length(Data.SeuObj@meta.data[["orig.ident"]]))
    }
    
    scRNA_SeuObj.list[[i]] <- Data.SeuObj
    names(scRNA_SeuObj.list)[[i]] <- list_files.df$Folder[i]
  }
  rm(i,j,Folder,Data.dgCMatrix,Data.SeuObj)

##### 01 Combine different datasets before QC  #####  

    # normalize and identify variable features for each dataset independently
  set.seed(1) # Fix the seed
  scRNA_SeuObj.list <- lapply(X = scRNA_SeuObj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  set.seed(1) # Fix the seed
  features <- SelectIntegrationFeatures(object.list = scRNA_SeuObj.list)
  
  ## Perform integration
  set.seed(1) # Fix the seed
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA_SeuObj.list, anchor.features = features)
  # this command creates an 'integrated' data assay
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- IntegrateData(anchorset = scRNA.anchors)
  
  set.seed(1) # Fix the seed
  DefaultAssay(scRNA.SeuObj) <- "integrated"
  
  save.image(paste0(Save.Path,"/01_Combine_different_datasets_before_QC.RData"))     
  
    
##### 02 Quality Control  #####
  dir.create(paste0(Save.Path,"/",ProjectName,"_QC"))
  ## QC for all samples
  scRNA.SeuObj_Ori <- scRNA.SeuObj # Save the original obj
  #Test# scRNA.SeuObj_Ori.list <- SplitObject(scRNA.SeuObj_Ori, split.by = "ID")
  scRNA.SeuObj_QCTry <- scRNAQC(scRNA.SeuObj,FileName = paste0(Version,"/",ProjectName,"_QC/",ProjectName,"_QCTry"))
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  
  rm(scRNA.anchors,scRNA.SeuObj)
  
  ## QC for each sample for the new integration
  scRNA_SeuObj_QC.list <- list()
  for (i in 1:length(scRNA_SeuObj.list)) {
    
    # scRNA_SeuObj.list[[i]] <- Data.SeuObj
    # names(scRNA_SeuObj_QC.list)[[i]] <- names
    Name <- names(scRNA_SeuObj.list)[[i]]
    
    scRNA_SeuObj_QC.list[[i]] <- scRNAQC(scRNA_SeuObj.list[[i]],
                                         FileName = paste0(Version,"/",ProjectName,"_QC/",ProjectName,"_", Name,"_QC"))
    names(scRNA_SeuObj_QC.list)[[i]] <- Name  

    }
  rm(i,Name)
  
  save.image(paste0(Save.Path,"/02_Quality_Control.RData"))      
  
  
##### 03 Combine different data sets after QC  #####

  # normalize and identify variable features for each dataset independently
  set.seed(1) # Fix the seed
  scRNA_SeuObj_QC.list <- lapply(X = scRNA_SeuObj_QC.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  set.seed(1) # Fix the seed
  features <- SelectIntegrationFeatures(object.list = scRNA_SeuObj_QC.list)
  
  ## Perform integration
  set.seed(1) # Fix the seed
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA_SeuObj_QC.list, anchor.features = features)
  # this command creates an 'integrated' data assay
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- IntegrateData(anchorset = scRNA.anchors)
  
  
  ## Check QC
  scRNAQC(scRNA.SeuObj,AddMitInf = "No",CheckOnly="Yes",FileName = paste0(Version,"/",ProjectName,"_QC/",ProjectName,"_QC_Check"))
  
  save.image(paste0(Save.Path,"/03_Combine_different_data_sets_after_QC.RData"))      
  
##### 04 Perform an integrated analysis #####
  # Run the standard workflow for visualization and clustering
  
  # # # !!
  # set.seed(1) # Fix the seed
  # all.genes <- rownames(scRNA.SeuObj)
  # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, features = all.genes)
  
  ## Issues: re_clustering in seurat v3
  ## https://github.com/satijalab/seurat/issues/1528
  # DefaultAssay(scRNA.SeuObj) <- "RNA"
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(scRNA.SeuObj) <- "integrated"
  
  ##?
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  
  # ## Run if use filter
  # set.seed(1) # Fix the seed
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  
  
  ### RunPCA
  # set.seed(1) # Fix the seed
  # scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 30, verbose = FALSE)
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, features = VariableFeatures(object = scRNA.SeuObj))
  
  print(scRNA.SeuObj[["pca"]], dims = 1:5, nfeatures = 5)
  
  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/PBMC_PCA.pdf"),
    width = 10,  height = 8
  )
  VizDimLoadings(scRNA.SeuObj, dims = 1:2, reduction = "pca")
  DimPlot(scRNA.SeuObj, reduction = "pca")
  DimHeatmap(scRNA.SeuObj, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(scRNA.SeuObj, dims = 1:15, cells = 500, balanced = TRUE)
  DimHeatmap(scRNA.SeuObj, dims = 16:30, cells = 500, balanced = TRUE)
  
  # # Determine the 'dimensionality' of the dataset
  # # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # # computation time
  # scRNA.SeuObj <- JackStraw(scRNA.SeuObj, num.replicate = 100)
  # scRNA.SeuObj <- ScoreJackStraw(scRNA.SeuObj, dims = 1:20)
  # JackStrawPlot(scRNA.SeuObj, dims = 1:20)
  ElbowPlot(scRNA.SeuObj, ndims = 50)
  dev.off()
  
  ElbowPlot(scRNA.SeuObj, ndims = 50)
  
  ## Issues: RunUMAP causes R exit
  ## https://github.com/satijalab/seurat/issues/2259
  # The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
  # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
  # This message will be shown once per session
  #### UMAP
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:30)
  
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)
  
  #### tSNE
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- RunTSNE(scRNA.SeuObj, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)
  
  
  ## Visualization
  DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
  
  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_Cluster.pdf"),
    width = 10,  height = 8
  )
  
  DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "sample") %>% 
    BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)
  
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE, label.size = 4) %>% 
    BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                   SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)
  
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>% 
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>% 
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
  
  ## tSNE
  DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = "sample") %>% 
    BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
  
  dev.off()
  # graphics.off()
  
  rm(PBMC.TN138_QC, PBMC.TN139_QC, PBMC.TN146_QC, PBMC.TN148_QC,scRNA.SeuObj_QCTry)
  
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))        
  
  ##### Meta Table  #####
  
  ## Before QC
  Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
  colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
  Meta.df[1,1] <- c("EO.M")  # TN138
  Meta.df[1,2] <- ncol(scRNA_SeuObj.list[[1]]@assays[["RNA"]]@counts)
  Meta.df[1,3] <- nrow(scRNA_SeuObj.list[[1]]@assays[["RNA"]]@counts)
  
  Meta.df[2,1] <- c("LO.M")  # TN139
  Meta.df[2,2] <- ncol(scRNA_SeuObj.list[[2]]@assays[["RNA"]]@counts)
  Meta.df[2,3] <- nrow(scRNA_SeuObj.list[[2]]@assays[["RNA"]]@counts)
  
  Meta.df[3,1] <- c("LO.F")  # TN146
  Meta.df[3,2] <- ncol(scRNA_SeuObj.list[[3]]@assays[["RNA"]]@counts)
  Meta.df[3,3] <- nrow(scRNA_SeuObj.list[[3]]@assays[["RNA"]]@counts)
  
  Meta.df[4,1] <- c("EO.F")  # TN148
  Meta.df[4,2] <- ncol(scRNA_SeuObj.list[[4]]@assays[["RNA"]]@counts)
  Meta.df[4,3] <- nrow(scRNA_SeuObj.list[[4]]@assays[["RNA"]]@counts)
  
  # Summary to Meta table
  Meta.df[5,1] <- c("Summary")
  Meta.df[5,2] <- ncol(scRNA.SeuObj_Ori@assays[["RNA"]]@counts)
  Meta.df[5,3] <- nrow(scRNA.SeuObj_Ori@assays[["RNA"]]@counts)
  
  ## After QC
  colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
  Meta.df[6,1] <- c("EO.M.QC")  # TN138
  Meta.df[6,2] <- ncol(scRNA_SeuObj_QC.list[[1]]@assays[["RNA"]]@counts)
  Meta.df[6,3] <- nrow(scRNA_SeuObj_QC.list[[1]]@assays[["RNA"]]@counts)
  
  Meta.df[7,1] <- c("LO.M.QC")  # TN139
  Meta.df[7,2] <- ncol(scRNA_SeuObj_QC.list[[2]]@assays[["RNA"]]@counts)
  Meta.df[7,3] <- nrow(scRNA_SeuObj_QC.list[[2]]@assays[["RNA"]]@counts)
  
  Meta.df[8,1] <- c("LO.F.QC")  # TN146
  Meta.df[8,2] <- ncol(scRNA_SeuObj_QC.list[[3]]@assays[["RNA"]]@counts)
  Meta.df[8,3] <- nrow(scRNA_SeuObj_QC.list[[3]]@assays[["RNA"]]@counts)
  
  Meta.df[9,1] <- c("EO.F.QC")  # TN148
  Meta.df[9,2] <- ncol(scRNA_SeuObj_QC.list[[4]]@assays[["RNA"]]@counts)
  Meta.df[9,3] <- nrow(scRNA_SeuObj_QC.list[[4]]@assays[["RNA"]]@counts)
  
  # Summary to Meta table
  Meta.df[10,1] <- c("Summary")
  Meta.df[10,2] <- ncol(scRNA.SeuObj@assays[["RNA"]]@counts)
  Meta.df[10,3] <- nrow(scRNA.SeuObj@assays[["RNA"]]@counts)
  
  
  write.table( Meta.df ,
               file = paste0(Save.Path,"/PBMC_CellCount_Meta.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )
  
  
##### 05 Identify conserved cell type markers  ##### 
  ## Identify conserved cell type markers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  set.seed(1) # Fix the seed
  PBMC.markers <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  library("magrittr")
  library("dplyr")
  
  # https://github.com/satijalab/seurat/issues/2960
  # Filter the top markers and plot the heatmap
  top_NSet = 7
  PBMC.markers %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  DoHeatmap(scRNA.SeuObj, features = top_N$gene) + NoLegend()
  write.table(top_N, file=paste0(Save.Path,"/PBMC_ClusterMarker_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  write.table(PBMC.markers, file=paste0(Save.Path,"/PBMC_ClusterMarker_AllGene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  
  pdf(
    file = paste0(Save.Path,"/PBMC_Heatmap_Cluster_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,size = 2,angle = 60) +
    scale_fill_gradient2(low="#5283ff",mid ="white", high ="#ff5c5c") +
    theme(axis.text.y = element_text(size  = 5)) +
    theme(legend.position = "bottom" )
  
  dev.off()
  
  
  # --------------- Check specific tissue marker --------------- #
  
  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_CTMarker.pdf"),
    width = 10,  height = 8
  )
  
  # PMID: 31771616 #!!!!!!
  FeaturePlot(scRNA.SeuObj, features = c("Cd3d", "Cd4", "Cd8a", "Csf1r", "Foxp3", "S100a9"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # T Cell: Cd3d;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r; regulatory T cells(Treg): Foxp3; Neutrophils: S100a9
  
  # PMID: 34296197 #!!!!!
  FeaturePlot(scRNA.SeuObj, features = c("Cd3d", "Cd3e", "Lyz1", "Lyz2","Clu","Cd79a","Ms4a1","Nkg7","Gzmb"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # T Cell: Cd3d,Cd3e;  Macrophages: Lyz; Mast Cell: Clu; B Cell: Cd79a,Ms4a1; NK Cell: Nkg7,Gzmb
  
  # http://biocc.hrbmu.edu.cn/CellMarker/
  # Mast cell
  FeaturePlot(scRNA.SeuObj, features = c("Cd117", "Cd25","Cd203c","Slc18a2","Kit","Fcer1a","Cd9"), min.cutoff = "q9", coord.fixed = 1)
  # PMID: 30356731
  FeaturePlot(scRNA.SeuObj, features = c("Cd9"), min.cutoff = "q9", coord.fixed = 1)
  # PMID: 34296197 #!!!!!
  FeaturePlot(scRNA.SeuObj, features = c("Cpa3"), min.cutoff = "q9", coord.fixed = 1)
  # https://www.panglaodb.se/markers.html?cell_type=%27Mast%20cells%27
  
  ## Mac
  # Macrophage-Markers
  # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
  ## M0
  FeaturePlot(scRNA.SeuObj, features = c("Cd68", "Adgre1","Cd14","Csf1r","Ly6c1",
                                          "Cx3cr1","Fcgr1a","Itgam","Mertk"), min.cutoff = "q9", coord.fixed = 1)
  
  ## M1
  # M1 http://biocc.hrbmu.edu.cn/CellMarker/
  FeaturePlot(scRNA.SeuObj, features = c("Cd16","Cd32","Cd64","Cd68","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)
  # M1 https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
  FeaturePlot(scRNA.SeuObj, features = c("Marco","Nos2","Tlr2","Cd80","Cd86","Csf2",
                                          "Tnf","Il1b","Il6","Tlr4","Cxcl2","Ifng","Il1r1"), min.cutoff = "q9", coord.fixed = 1)
  FeaturePlot(scRNA.SeuObj, features = c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)
  
  
  ## M2
  # M2 http://biocc.hrbmu.edu.cn/CellMarker/
  FeaturePlot(scRNA.SeuObj, features = c("Chil3","Csf1r","Mrc1","Pparg","Arg1","Cd163","Clec10a","Clec7a",
                                          "Cd206","Cd209","Ccl18","Fizz1"), min.cutoff = "q9", coord.fixed = 1)
  # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
  FeaturePlot(scRNA.SeuObj, features = c("Cd115", "Cd206", "Pparg", "Arg1", "Cd163", "Cd301", 
                                          "Dectin-1", "Pdcd1lg2", "Fizz1"), min.cutoff = "q9", coord.fixed = 1)
  
  FeaturePlot(scRNA.SeuObj, features = c("Chil3"), min.cutoff = "q9", coord.fixed = 1)
  
  
  
  ## Tumor associated macrophage(TAM)
  FeaturePlot(scRNA.SeuObj, features = c("Ccr2","Csf1r","Marco","Pdl2","Cd40","Ccl2","Csf1","Cd16"), 
              min.cutoff = "q9", coord.fixed = 1)
  
  # Erythrocytes
  FeaturePlot(scRNA.SeuObj, features = c("Hbb-bs"), min.cutoff = "q9", coord.fixed = 1)
  # Platelet 
  FeaturePlot(scRNA.SeuObj, features = c("Ppbp"), min.cutoff = "q9", coord.fixed = 1)
  
  ## Summary
  markers.to.plot <- c("Cd3d","Cd3e", "Cd4","Cd8a", "Csf1r", "Lyz2","Chil3","Il1b", "S100a9","Nkg7",
                       "Gzmb", "Cd79a", "Ms4a1","Clu","Hbb-bs","Ppbp")
  
  FeaturePlot(scRNA.SeuObj, features = markers.to.plot, min.cutoff = "q9", coord.fixed = 1)
  # T Cell: Cd3d,Cd3e;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r,Lyz,Chil3;  Neutrophils: S100a9; 
  # NK Cell: Nkg7,Gzmb; B Cell: Cd79a,Ms4a1; Mast Cell: Clu; Erythrocytes: Hbb-bs; Platelet: Ppbp
  
  
  dev.off()
  
  save.image(paste0(Save.Path,"/05_Identify_conserved_cell_type_markers.RData"))  
  
##### 06 Cell type annotation  #####
  # scRNA.SeuObj.copy <- scRNA.SeuObj
  
  ## CD4+T: CD4+T Cell; CD8+T: CD8+T Cell; T: T Cell; B: B Cell; Mac: Macrophages;
  ## Neu: Neutrophils; NK: NK Cell; Mast: Mast Cell; Ery: Erythrocytes;
  ## Thr: Thrombocytes
  scRNA.SeuObj <- RenameIdents(scRNA.SeuObj, `0` = "CD4+T", `1` = "B", `2` = "Mac3",
                                `3` = "Neu", `4` = "CD8+T", `5` = "CD8+T", `6` = "Mac2", `7` = "CD4+T", 
                                `8` = "NK", `9` = "Neu",`10` = "Mast1", `11` = "T", `12` = "Ery", `13` = "Mac1", 
                                `14` = "B", `15` = "B", `16` = "Mast2", `17` = "Mac0", `18` = "Neu")
  
  Cell_Type_Order.set <- c("T", "CD4+T", "CD8+T", "B" , "Mac0", "Mac1", "Mac2", "Mac3", "Mast1", "Mast2", "NK", "Neu", "Ery")
  
  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  # Idents(scRNA.SeuObj) <- "celltype"
  
  ## Heatmap
  Heatmap_Color.lt <- list(low="#5283ff",mid ="white", high ="#ff5c5c")
  
  pdf(
    file = paste0(Save.Path,"/PBMC_Heatmap_CellType_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,size = 3,angle = 60) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]], 
                         high = Heatmap_Color.lt[["high"]]) +
    theme(axis.text.y = element_text(size  = 5)) +
    theme(legend.position = "bottom")+
    theme(aspect.ratio=1)
  dev.off()
  
  # ## Ch
  # # https://github.com/satijalab/seurat/issues/2960
  # set.seed(1) # Fix the seed
  # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,size = 4,angle = 90) + NoLegend()
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,group.by = "celltype",size = 4,angle = 90) + NoLegend()
  
  # Color
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,group.by = "celltype",size = 3,angle = 60) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]], 
                         high = Heatmap_Color.lt[["high"]]) 
  
  
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,group.by = "seurat_clusters",size = 3,angle = 90) + NoLegend()
  
  pdf( 
    file = paste0(Save.Path,"/PBMC_Heatmap_CellType_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(scRNA.SeuObj, features = top_N$gene,group.by = "celltype",size = 2,angle = 45) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]], 
                         high = Heatmap_Color.lt[["high"]])  +
    theme(axis.text.y = element_text(size = 5)) +
    theme(legend.position = "bottom" )
  
  dev.off()
  
  ## UMAP tSNE
  DimPlot(scRNA.SeuObj, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))
  
  DimPlot(scRNA.SeuObj,group.by = "celltype",label.size = 7,label = TRUE,  
          pt.size =2) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_CellType")
  DimPlot(scRNA.SeuObj,group.by = "sample",  
          pt.size =0.5) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_Sample")
  DimPlot(scRNA.SeuObj,group.by = "seurat_clusters",label.size = 7, label = TRUE,  
          pt.size =1) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_Clusters")
  
  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_CellType_Sup.pdf"),
    width = 10,  height = 8
  )
  ##
  
  
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE) %>% 
    BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                   SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)
  
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>% 
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
  
  DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>% 
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
  
  ## tSNE
  DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = "sample") %>% 
    BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
  
  dev.off()
  
  
  ## DotPlot
  DotPlot_Color1.set <- c("#de3767", "#de3767", "#4169e1", "#4169e1")
  DotPlot_Color2.set <- c("#5b8e7d","#7b2cbf")
  DotPlot_Color3.set <- c("#de3767", "#4169e1")
  
  pdf( 
    file = paste0(Save.Path,"/PBMC_DotPlot_CellType",".pdf"),
    width = 10,  height = 8
  )
  
  # https://satijalab.org/seurat/reference/dotplot
  DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = c("lightgrey", "blue"), 
          dot.scale = 8) + RotatedAxis()%>% 
    BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1, TitleSize = 20, xangle =90,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1,XtextSize=12,  YtextSize=12)
  
  DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = DotPlot_Color1.set, 
          dot.scale = 8, split.by = "sample") + RotatedAxis()
  
  
  # https://github.com/satijalab/seurat/issues/1541
  DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = DotPlot_Color2.set, 
          dot.scale = 8, split.by = "Cachexia") + RotatedAxis()
  
  DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = DotPlot_Color3.set, 
          dot.scale = 8, split.by = "Sex") + RotatedAxis()
  
  dev.off()
  
  rm(top_N, top_NSet)
  
  save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))
  
  # ##### Export marker gene from specific cluster #####
  #   # For performing differential expression after integration, we switch back to the original data
  #   set.seed(1) # Fix the seed
  #   DefaultAssay(scRNA.SeuObj) <- "RNA"
  #   
  #   # nk.markers <- FindConservedMarkers(scRNA.SeuObj, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
  #   library(BiocManager)
  #   library(multtest)
  #   nk.markers <- FindConservedMarkers(scRNA.SeuObj, ident.1 = 'NK', grouping.var = "sample", verbose = FALSE)
  #   head(nk.markers)
  #   
  #   rm(nk.markers)
  
  # ##### Identify differential expressed genes across conditions  #####
  #   library(ggplot2)
  #   library(cowplot)
  #   theme_set(theme_cowplot())
  #   CD4T.cells <- subset(scRNA.SeuObj, idents = "CD4+T")
  #   Idents(CD4T.cells) <- "Cachexia"
  #   avg.CD4T.cells <- as.data.frame(log1p(AverageExpression(CD4T.cells, verbose = FALSE)$RNA))
  #   avg.CD4T.cells$gene <- rownames(avg.CD4T.cells)
  #   
  #   MacrophageM2 <- subset(scRNA.SeuObj, idents = "Mac2")
  #   Idents(MacrophageM2) <- "Cachexia"
  #   avg.MacrophageM2 <- as.data.frame(log1p(AverageExpression(MacrophageM2, verbose = FALSE)$RNA))
  #   avg.MacrophageM2$gene <- rownames(avg.MacrophageM2)
  #   
  #   genes.to.label = c("Sox17", "Mrpl15", "Lypla1", "Tcea1", "Rgs20", "Atp6v1h", "Rb1cc1", "4732440D04Rik", "St18")
  #   p1 <- ggplot(avg.CD4T.cells, aes(EO, LO)) + geom_point() + ggtitle("Cachexia T.cells")
  #   p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  #   p2 <- ggplot(avg.MacrophageM2, aes(EO, LO)) + geom_point() + ggtitle("Cachexia Macrophage")
  #   p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
  #   p1 + p2
  #   rm(p1 , p2 ,CD4T.cells, MacrophageM2, avg.CD4T.cells, avg.MacrophageM2)
  
  
##### 07 Count Cell number  #####
  Pheno.df <- data.frame(sample = scRNA.SeuObj@meta.data[["sample"]],celltype = scRNA.SeuObj@meta.data[["celltype"]],
                         Cachexia = scRNA.SeuObj@meta.data[["Cachexia"]],Sex = scRNA.SeuObj@meta.data[["Sex"]])
  # Pheno.df.table <- table(Pheno.df) %>% as.data.frame()
  
  Freq_sample.df <- table(Pheno.df$sample) %>% as.data.frame()
  Freq_CT.df <- table(Pheno.df$celltype) %>% as.data.frame()
  Freq_Cach.df <- table(Pheno.df$Cachexia) %>% as.data.frame()
  Freq_Sex.df <- table(Pheno.df$Sex) %>% as.data.frame()
  
  ##
  Pheno_EO_M.df <- Pheno.df[Pheno.df$sample=="EO.M",]
  Pheno_LO_M.df <- Pheno.df[Pheno.df$sample=="LO.M",]
  Pheno_EO_F.df <- Pheno.df[Pheno.df$sample=="EO.F",]
  Pheno_LO_F.df <- Pheno.df[Pheno.df$sample=="LO.F",]
  
  # Count EO_M_CT
  Freq_EO_M_CT.df <- table(Pheno_EO_M.df$celltype) %>% as.data.frame()
  Freq_EO_M_CT.df <- data.frame(Type="EO.M",Freq_EO_M_CT.df)
  Freq_EO_M_CT.df$Percent <- Freq_EO_M_CT.df$Freq/sum(Freq_EO_M_CT.df$Freq)
  
  # Count LO_M_CT
  Freq_LO_M_CT.df <- table(Pheno_LO_M.df$celltype) %>% as.data.frame()
  Freq_LO_M_CT.df <- data.frame(Type="LO.M",Freq_LO_M_CT.df)
  Freq_LO_M_CT.df$Percent <- Freq_LO_M_CT.df$Freq/sum(Freq_LO_M_CT.df$Freq)
  
  # Count EO_F_CT
  Freq_EO_F_CT.df <- table(Pheno_EO_F.df$celltype) %>% as.data.frame()
  Freq_EO_F_CT.df <- data.frame(Type="EO.F",Freq_EO_F_CT.df)
  Freq_EO_F_CT.df$Percent <- Freq_EO_F_CT.df$Freq/sum(Freq_EO_F_CT.df$Freq)
  
  # Count LO_F_CT
  Freq_LO_F_CT.df <- table(Pheno_LO_F.df$celltype) %>% as.data.frame()
  Freq_LO_F_CT.df <- data.frame(Type="LO.F",Freq_LO_F_CT.df)
  Freq_LO_F_CT.df$Percent <- Freq_LO_F_CT.df$Freq/sum(Freq_LO_F_CT.df$Freq)
  
  # Combind all count of sample
  Freq_All.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                       Freq_EO_F_CT.df,Freq_LO_F_CT.df)
  Freq_All.df <- data.frame(Index = row.names(Freq_All.df),Freq_All.df )
  colnames(Freq_All.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")
  # Freq_All.df$Index <- factor(Freq_All.df$Index,
  #                                              levels = Freq_All.df$Index)
  
  #### LinePlot ####  
  # https://ithelp.ithome.com.tw/articles/10186047
  # Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
  #                                    levels = sort(unique(as.character(Freq_All.df$Cell_Type))))
  
  Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
                                  levels = Cell_Type_Order.set)
  
  CellNum_P1 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Number, 
                                        colour = Pheno_Type, group = Pheno_Type)) + 
    geom_line(linetype = "dashed",size=1.5) + 
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  CellNum_P1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=15,  YtextSize=15, xangle = 90,
                                LegTextSize = 15) + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P1
  CellNum_P1
  
  CellNum_P2 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Percent, 
                                        colour = Pheno_Type, group = Pheno_Type)) + 
    geom_line(linetype = "dashed",size=1.5) + 
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
  
  CellNum_P2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18, xangle = 90,
                                LegTextSize = 15) + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P2
  CellNum_P2
  
  
  #### All type compare to Combine Sex ####  
  ##
  Pheno_EO.df <- Pheno.df[Pheno.df$Cachexia=="EO",]
  Pheno_LO.df <- Pheno.df[Pheno.df$Cachexia=="LO",]
  
  # Count EO_CT
  Freq_EO_CT.df <- table(Pheno_EO.df$celltype) %>% as.data.frame()
  Freq_EO_CT.df <- data.frame(Type="EO",Freq_EO_CT.df)
  Freq_EO_CT.df$Percent <- Freq_EO_CT.df$Freq/sum(Freq_EO_CT.df$Freq)
  
  # Count LO_CT
  Freq_LO_CT.df <- table(Pheno_LO.df$celltype) %>% as.data.frame()
  Freq_LO_CT.df <- data.frame(Type="LO",Freq_LO_CT.df)
  Freq_LO_CT.df$Percent <- Freq_LO_CT.df$Freq/sum(Freq_LO_CT.df$Freq)
  
  # Combind all count of sample
  Freq_All_Ca.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                          Freq_EO_F_CT.df,Freq_LO_F_CT.df,
                          Freq_EO_CT.df,Freq_LO_CT.df)
  
  Freq_All_Ca.df <- data.frame(Index = row.names(Freq_All_Ca.df),Freq_All_Ca.df )
  colnames(Freq_All_Ca.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")
  
  # Change the order
  # https://blog.csdn.net/weixin_48172266/article/details/117537465
  # CTOrder.set <- factor(Freq_All_Ca.df$Cell_Type,
  #                       levels = sort(unique(as.character(Freq_All_Ca.df$Cell_Type))))
  # 
  # Freq_All_Ca.df <- Freq_All_Ca.df %>% 
  #                   mutate(Cell_Type = CTOrder.set)
  
  ## Freq_All_Ca.df$Cell_Type <- factor(Freq_All_Ca.df$Cell_Type,
  ##                                    levels = sort(unique(as.character(Freq_All_Ca.df$Cell_Type))))
  
  # Freq_All_Ca.df <- read.csv(paste0(Save.Path,"/PBMC_CT_Count_CH.txt"),sep = "\t")
  Freq_All_Ca.df$Percent <- as.numeric(Freq_All_Ca.df$Percent)
  
  write.table( Freq_All_Ca.df ,
               file = paste0(Save.Path,"/PBMC_CellCount_CT_Ca.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )
  
  # https://stackoverflow.com/questions/27350243/ggplot-line-graph-with-different-line-styles-and-markers/27350366
  # https://www.coder.work/article/6971741
  # https://stackoverflow.com/questions/11344561/controlling-line-color-and-line-type-in-ggplot-legend
  
  #### LinePlot ####   
  # Freq_All_Ca.df$Cell_Type <- factor(Freq_All_Ca.df$Cell_Type,
  #                                    levels = sort(unique(as.character(Freq_All_Ca.df$Cell_Type))))
  Freq_All_Ca.df$Cell_Type <- factor(Freq_All_Ca.df$Cell_Type,
                                     levels = Cell_Type_Order.set)
  
  CellNum_P3 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Number, 
                                           colour = Pheno_Type,
                                           group = Pheno_Type,linetype=Pheno_Type
  )) + 
    geom_line(size=1.5) + 
    scale_linetype_manual(name="Pheno_Type", 
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) + 
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") + theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  CellNum_P3 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.4, 0.8),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18,xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P3
  CellNum_P3
  
  library(eoffice)
  topptx(CellNum_P3,paste0(Save.Path,"/Temp.pptx"))
  
  CellNum_P4 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Percent, 
                                           colour = Pheno_Type,
                                           group = Pheno_Type,linetype=Pheno_Type
  )) + 
    geom_line(size=1.5) + 
    scale_linetype_manual(name="Pheno_Type", 
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) + 
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
    #theme_set(theme_bw())+ # Remove the background
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  CellNum_P4 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.15, 0.82),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P4
  CellNum_P4
  
  ##### Clean the dataframe #####    
  CellPheno.lt <- list(Pheno.df, Pheno_EO.df, Pheno_LO.df,
                       Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)
  names(CellPheno.lt) <- c("Pheno.df", "Pheno_EO.df", "Pheno_LO.df",
                           "Pheno_EO_M.df","Pheno_LO_M.df","Pheno_EO_F.df","Pheno_LO_F.df")
  rm(Pheno_EO.df, Pheno_LO.df, Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)
  
  CellFreq.lt <- list(Freq_All_Ca.df, Freq_All.df, Freq_sample.df, 
                      Freq_Sex.df,Freq_Cach.df, Freq_CT.df, 
                      Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
                      Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)
  names(CellFreq.lt) <- c("Freq_All_Ca.df", "Freq_All.df", "Freq_sample.df", 
                          "Freq_Sex.df", "Freq_Cach.df", "Freq_CT.df", 
                          "Freq_EO_CT.df", "Freq_EO_F_CT.df", "Freq_EO_M_CT.df",
                          "Freq_LO_CT.df", "Freq_LO_F_CT.df", "Freq_LO_M_CT.df")
  rm(Freq_sample.df, 
     Freq_Sex.df,Freq_Cach.df, Freq_CT.df, 
     Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
     Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)
  
  ##### BarPlot #####  
  # https://blog.gtwang.org/r/ggplot2-tutorial-layer-by-layer-plotting/3/
  colnames(Pheno.df) <- c("sample","Cell_Type","Cachexia","Sex")
  Pheno.df$Cell_Type <- factor(Pheno.df$Cell_Type,
                               levels = sort(unique(as.character(Pheno.df$Cell_Type))))
  
  # sample
  BarPlot1_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) + 
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_1
  BarPlot1_1
  
  BarPlot1_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) + 
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot1_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_2
  BarPlot1_2
  
  # Cachexia
  BarPlot2_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) + 
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot2_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_1
  BarPlot2_1
  
  BarPlot2_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) + 
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot2_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_2
  BarPlot2_2
  
  # Sex
  BarPlot3_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) + 
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot3_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_1
  BarPlot3_1
  BarPlot3_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) + 
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot3_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  + 
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_2
  
  BarPlot3_2
  
  
  ##### Export PDF file #####
  pdf(file = paste0(Save.Path,"/PBMC_CellCount_LinePlot.pdf"),
      width = 7, height = 7 )
  CellNum_P4
  CellNum_P3
  CellNum_P1
  CellNum_P2
  BarPlot1_1
  BarPlot1_2
  BarPlot2_1
  BarPlot2_2
  BarPlot3_1
  BarPlot3_2
  dev.off() # graphics.off()
  
  rm(CellNum_P1, CellNum_P2, CellNum_P3, CellNum_P4, BarPlot1_1, BarPlot1_2,
     BarPlot2_1, BarPlot2_2, BarPlot3_1, BarPlot3_2)
  
  save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))
  
  
##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  scRNA.SeuObj$celltype.Cachexia <- paste(Idents(scRNA.SeuObj), scRNA.SeuObj$Cachexia, sep = "_")
  scRNA.SeuObj$celltype.Cachexia.gender <- paste(Idents(scRNA.SeuObj), scRNA.SeuObj$Cachexia, scRNA.SeuObj$Sex, sep = "_")
  Idents(scRNA.SeuObj) <- "celltype.Cachexia.gender"
  
  scRNA.SeuObj$Cachexia.gender <- paste(scRNA.SeuObj$Cachexia, scRNA.SeuObj$Sex, sep = "_")
  
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  
  ##-------------- Find Marker gene in Male --------------##
  CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  # CellType.list <- CellType.list[-9] 
  
  dir.create(paste0(Save.Path,"/PBMC_SSA_Male_FindMarkers"))
  
  # About 15 mins
  CCMarker_Male.lt <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      CCMarker_Male.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                            paste0(CellType.list[i],"_EO_Male"), 
                                            paste0(CellType.list[i],"_LO_Male"),
                                            CellType.list[i],
                                            Path = Save.Path,
                                            ResultFolder = "PBMC_SSA_Male_FindMarkers")
      # names(CCMarker_Male.lt)[[i]] <- paste0("CCMarker_Male.lt.",CellType.list[i])
      names(CCMarker_Male.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)
  
  CCMarker_Male.lt <- CCMarker_Male.lt[!unlist(lapply(CCMarker_Male.lt,is.null))]
  
  
  ## Generate pdf and tif file for Male VolcanoPlot 
  dir.create(paste0(Save.Path,"/PBMC_SSA_Male_VolcanoPlot/"))
  
  pdf(file = paste0(Save.Path,"/PBMC_SSA_Male_VolcanoPlot/PBMC_SSA_Male_VolcanoPlot.pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+ 
              ggtitle(paste0("PBMC_Male_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/PBMC_SSA_Male_VolcanoPlot/",CellType.list[i],".tif"), 
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]])+ 
              ggtitle(paste0("PBMC_Male_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i)
  
  ##-------------- Find Marker gene in Female --------------##
  CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  # CellType.list <- CellType.list[-9] 
  
  dir.create(paste0(Save.Path,"/PBMC_SSA_Female_FindMarkers"))
  
  # About 15 mins
  CCMarker_Female.lt <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      CCMarker_Female.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                              paste0(CellType.list[i],"_EO_Female"), 
                                              paste0(CellType.list[i],"_LO_Female"),
                                              CellType.list[i],
                                              Path = Save.Path,
                                              ResultFolder = "PBMC_SSA_Female_FindMarkers")
      # names(CCMarker_Female.lt)[[i]] <- paste0("CCMarker_Female.lt.",CellType.list[i])
      names(CCMarker_Female.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)
  
  CCMarker_Female.lt <- CCMarker_Female.lt[!unlist(lapply(CCMarker_Female.lt,is.null))]
  
  ## Generate pdf and tif file for Female VolcanoPlot 
  dir.create(paste0(Save.Path,"/PBMC_SSA_Female_VolcanoPlot/"))
  
  pdf(file = paste0(Save.Path,"/PBMC_SSA_Female_VolcanoPlot/PBMC_SSA_Female_VolcanoPlot.pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+ 
              ggtitle(paste0("PBMC_Female_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/PBMC_SSA_Female_VolcanoPlot/",CellType.list[i],".tif"), 
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]])+ 
              ggtitle(paste0("PBMC_Female_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i)
  
  save.image(paste0(Save.Path,"/08_1_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SSA).RData"))
  
  
##### 08_1 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ##-------------- Intersect_CellType --------------##
  #CCMarker_Male_Ori.lt <- CCMarker_Male.lt
  #CCMarker_Female_Ori.lt <- CCMarker_Female.lt   
  #CellType_Ori.list <- CellType.list
  
  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))
  
  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]
  
  CellType.list <- names(CCMarker_Male.lt)
  
  ##-------------- Venn Pos --------------##
  source("FUN_Venn.R")
  # pdf(file = paste0(Save.Path,"/PBMC_Female_VolcanoPlot.pdf"),width = 7, height = 7 )
  
  dir.create(paste0(Save.Path,"/PBMC_SSA_VennDiagrame"))
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                               CellType.list[i],"Pos","#9d0208","#f08080",SampleType="PBMC",
                                               PathName = paste0(Save.Path,"/PBMC_SSA_VennDiagrame"))
      names(Venn_CCMarker_Pos)[[i]] <- paste0("Venn_CCMarker.",CellType.list[i],"_Pos")
    })
  }
  rm(i)
  
  ##-------------- Venn Neg --------------##
  Venn_CCMarker_Neg <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      Venn_CCMarker_Neg[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.S_Neg_List"]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.S_Neg_List"]],
                                               CellType.list[i],"Neg","#00296b","#1368aa",SampleType="PBMC",
                                               PathName = paste0(Save.Path,"/PBMC_SSA_VennDiagrame"))
      
      names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_CCMarker.",CellType.list[i],"_Neg")
    })
  }
  rm(i)
  
  save.image(paste0(Save.Path,"/08_1_Find_CCmarker_in_different_Cell_type_and_VennDiagrame(SSA_IntersectCT).RData"))
  
  
##### 08_2 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  
  Idents(scRNA.SeuObj) <- "celltype.Cachexia"
  #CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  dir.create(paste0(Save.Path,"/PBMC_SPA_FindMarkers"))
  
  CCMarker_SPA.lt <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      CCMarker_SPA.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                           paste0(CellType.list[i],"_EO"), 
                                           paste0(CellType.list[i],"_LO"),
                                           CellType.list[i],
                                           Path = Save.Path,
                                           ResultFolder = "PBMC_SPA_FindMarkers")
      
      # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
      names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)
  
  CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]
  
  
  ## Generate pdf and tif file for VolcanoPlot 
  dir.create(paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/"))
  
  pdf(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot.pdf"),width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+ 
              ggtitle(paste0("PBMC_",CellType.list[i]))
      )
    })
  }
  # graphics.off()
  dev.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/PBMC_SPA_VolcanoPlot/PBMC_SPA_VolcanoPlot",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]])+ ggtitle(paste0("PBMC_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i)
  
  save.image(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))
  
  
  #####------------------------------------------------------------------------------------------------------------#####
  
##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_ggplot.R")
  
  # Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
                             col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
                             header = F,sep = "\t")
  
  # Convert Human gene to mouse
  Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
  for (i in 1:nrow(Pathway.all)) {
    #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t() 
    colnames(PathwayN)="Test"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
    Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
    
  }
  
  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
  
  save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))
  
##### 09_1 GSEA Analysis (SPA) #####
  GSEA_Large <- list()
  GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large.df.TOP <- GSEA_Large.df
  
  dir.create(paste0(Save.Path,"/PBMC_GSEA"))
  
  
  pdf(file = paste0(Save.Path, "/PBMC_GSEA/PBMC_GSEA_SPA.pdf"),width = 15, height = 7 )
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    
    fgseaRes <- GSEA_Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA_Large.Output[["Pathway.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA_Large.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large[[i]] <- Sum
    names(GSEA_Large)[[i]] <- paste0(CellType.list[i])
    
    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  
  dev.off()
  
  ## GSEA_Large.Sum.TOP ##
  GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
  GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Pathway_LargeTOP_SPA.txt"),sep="\t",
              row.names=F, quote = FALSE)
  
  ##### Bubble plot #####
  library(ggplot2)
  library(scales)
  GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")
  
  GSEA_Large.Sum.TOP$PhenoType <- factor(GSEA_Large.Sum.TOP$PhenoType,
                                         levels = Cell_Type_Order.set)
  
  GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
  GSEA_Large.Sum.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
  # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]
  
  pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SPA.lt[["BBPlot"]]
  GSEA_ggplot_SPA.lt[["BBPlot2"]]
  GSEA_ggplot_SPA.lt[["BBPlotB1"]]
  GSEA_ggplot_SPA.lt[["BBPlotB1"]]
  dev.off()
  
  
  ##### Extract SubType #####
  ## T Cell
  # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
  GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]
  
  BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]], 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  BBPlot_T
  
  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                           XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  # BBPlot_TB <- BBPlot_TB +theme(axis.title.y=element_blank(),
  #                  axis.text.y=element_blank(),
  #                  axis.ticks.y=element_blank())
  BBPlot_TB
  
  BBPlot_TB1 <- BBPlot_TB %>%
    insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_TB1
  
  
  pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB
  BBPlot_TB1
  dev.off()
  
  
  ## Mac
  GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]
  
  BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_size_area(max_size = 5)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]], 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  BBPlot_Mac
  
  BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                               XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  
  BBPlot_MacB1 <- BBPlot_MacB %>%
    insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_MacB1
  
  pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA_SubType_Mac.pdf"),width = 17, height = 20 )
  BBPlot_MacB
  BBPlot_MacB1
  dev.off()
  
  rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
     df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)
  
  ##### save.image #####
  save.image(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))    
  
##### 09_2 GSEA Analysis (SSA_MAle) #####
  GSEA_Large_Male <- list()
  GSEA_Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df
  
  
  pdf(file = paste0(Save.Path, "/PBMC_GSEA/PBMC_GSEA_SSA_Male.pdf"),width = 15, height = 7 )
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    
    #GSEA_Large_Male.Output <- FUN_GSEA_Large_MaleGeneSet(ranks,Pathway.all,10)
    GSEA_Large_Male.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    
    fgseaRes <- GSEA_Large_Male.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA_Large_Male.Output[["Pathway.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA_Large_Male.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large_Male[[i]] <- Sum
    names(GSEA_Large_Male)[[i]] <- paste0(CellType.list[i])
    
    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df <- rbind(GSEA_Large_Male.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df.TOP <- rbind(GSEA_Large_Male.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  
  dev.off()
  
  ## GSEA_Large_Male.Sum.TOP ##
  GSEA_Large_Male.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP)
  GSEA_Large_Male.Sum.TOP <- GSEA_Large_Male.Sum.TOP[,!colnames(GSEA_Large_Male.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large_Male.Sum.TOP, file=paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
              row.names=F, quote = FALSE)
  
  
  
  ##### Bubble plot #####
  library(ggplot2)
  library(scales)
  
  GSEA_Large_Male.Sum.TOP$PhenoType <- factor(GSEA_Large_Male.Sum.TOP$PhenoType,
                                              levels = Cell_Type_Order.set)
  
  GSEA_ggplot_SSA_Male.lt <- GSEA_ggplot(GSEA_Large_Male.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
  GSEA_Large_Male.Sum.TOP.S <- GSEA_ggplot_SSA_Male.lt[["GSEA_TOP.df"]]
  
  pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SSA_Male.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SSA_Male.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot2"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
  dev.off()
  
  
  ##### Extract SubType #####
  
  ## T Cell
  # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
  GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]
  
  BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]], 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  BBPlot_T
  
  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                           XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  BBPlot_TB
  
  BBPlot_TB1 <- BBPlot_TB %>%
    insert_left(GSEA_ggplot_SSA_Male.lt[["Y_Order"]],width = 0.2)
  BBPlot_TB1
  
  
  pdf(file = paste0(Save.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SSA_Male_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB
  BBPlot_TB1
  dev.off()
  
  
  
  
  
  