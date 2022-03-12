##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

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
  source("FUN_ReadscRNA.R")
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

##### Current path and new folder setting*  ##### 
  ProjectName = "CC"
  Sampletype = "PBMC"
  Version = paste0(Sys.Date(),"_","CC_PBMC")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)
  
  ## Import information
  InputFolder = "Input_files_10x" 
  InputAnno = "PBMC_Ano.csv"


##### Parameter setting* #####
  ClassSet1 = 1
  ClassSet2 = "Cachexia"
  ClassSet3 = "Sex"
  
##### Load datasets  #####
  ## Annotation table
  list_files.df <- read.csv(paste0(InputFolder,"/",InputAnno))
  Feature.set <- colnames(list_files.df)[-1]
  
  ## Read 10x files
  scRNA_SeuObj.list <- ReadscRNA(InputFolder, list_files.df, Mode="10x")
  
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
  ## Creative QC folder
  dir.create(paste0(Save.Path,"/",ProjectName,"_QC"))
  
  ## QC for all samples
  scRNA_Ori.SeuObj <- scRNA.SeuObj # Save the original obj
  #Test# scRNA.SeuObj_Ori.list <- SplitObject(scRNA.SeuObj_Ori, split.by = "ID")
  scRNA.SeuObj_QCTry <- scRNAQC(scRNA.SeuObj,FileName = paste0(Version,"/",ProjectName,"_QC/",ProjectName,"_QCTry"))
  rm(scRNA.anchors,scRNA.SeuObj,scRNA.SeuObj_QCTry)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay

  ## QC for each sample for the new integration
  scRNA_SeuObj_QC.list <- list()
  for (i in 1:length(scRNA_SeuObj.list)) {
    
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
  ## Creative Clusters folder
  dir.create(paste0(Save.Path,"/",ProjectName,"_Clusters"))
  
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
    file = paste0(Save.Path,"/",ProjectName,"_Clusters/",ProjectName,"_PCA.pdf"),
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
  DimPlot(scRNA.SeuObj, reduction = "umap", group.by = colnames(list_files.df)[3] ) %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
  
  pdf(
    file = paste0(Save.Path,"/",ProjectName,"_Clusters/",ProjectName,"_nlDR_Cluster.pdf"),
    width = 12,  height = 8
  )
    
    DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>% 
      BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)
  
   
    for (i in 1:(ncol(list_files.df)-1)) {
      print(DimPlot(scRNA.SeuObj, reduction = "umap", group.by = colnames(list_files.df)[i+1]) %>% 
            BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+ 
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2, split.by = colnames(list_files.df)[i+1], label = TRUE, label.size = 4) %>% 
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
      
      print(DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = colnames(list_files.df)[i+1]) %>% 
              BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+ 
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "tsne", ncol = 2, split.by = colnames(list_files.df)[i+1], label = TRUE, label.size = 4) %>% 
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
      
    }
    rm(i)

  dev.off()
  # graphics.off()
  
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))        
  
  ##### Meta Table  #####

    Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(Meta.df) <- c("Folder","Cell_Num","Gene_Num")

    ## Before QC
    for (i in 1:length(scRNA_SeuObj.list)) {
      Meta.df[i,1] <- names(scRNA_SeuObj.list[i])
      Meta.df[i,2] <- ncol(scRNA_SeuObj.list[[i]]@assays[["RNA"]]@counts)
      Meta.df[i,3] <- nrow(scRNA_SeuObj.list[[i]]@assays[["RNA"]]@counts)
      
    }

      # Summary to Meta table
      Meta.df[i+1,1] <- c("Summary")
      Meta.df[i+1,2] <- ncol(scRNA.SeuObj_Ori@assays[["RNA"]]@counts)
      Meta.df[i+1,3] <- nrow(scRNA.SeuObj_Ori@assays[["RNA"]]@counts)
    
    ## After QC
    for (j in 1:length(scRNA_SeuObj_QC.list)) {
      Meta.df[i+j+1,1] <- paste0(names(scRNA_SeuObj_QC.list[j]),".QC")
      Meta.df[i+j+1,2] <- ncol(scRNA_SeuObj_QC.list[[j]]@assays[["RNA"]]@counts)
      Meta.df[i+j+1,3] <- nrow(scRNA_SeuObj_QC.list[[j]]@assays[["RNA"]]@counts)
      
    }
    
      # Summary to Meta table
      Meta.df[i+j+2,1] <- c("Summary.QC")
      Meta.df[i+j+2,2] <- ncol(scRNA.SeuObj@assays[["RNA"]]@counts)
      Meta.df[i+j+2,3] <- nrow(scRNA.SeuObj@assays[["RNA"]]@counts)
      
    
    rm(i,j)
    
    Meta.df <- left_join(Meta.df,list_files.df, by = "Folder")
    Meta.df[is.na(Meta.df)] <- ""

  ## Export data        
  write.table( Meta.df ,
               file = paste0(Save.Path,"/",ProjectName,"_CellCount_Meta.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )
  


################## (Pending) Cell Cycle Regression ##################    
################## (Pending) Auto Cell type annotation ##################     
##### 05 Identify conserved cell type markers  ##### 
  ## Creative Cell type folder
  dir.create(paste0(Save.Path,"/",ProjectName,"_CT"))
  
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
  

  write.table(top_N, file=paste0(Save.Path,"/",ProjectName,"_Clusters/",ProjectName,"_ClusterMarker_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  write.table(PBMC.markers, file=paste0(Save.Path,"/",ProjectName,"_Clusters/",ProjectName,"_ClusterMarker_AllGene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  
  pdf(
    file = paste0(Save.Path,"/",ProjectName,"_Clusters/",ProjectName,"_Heatmap_Cluster_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
    DoHeatmap(scRNA.SeuObj, features = top_N$gene,size = 2,angle = 60) +
      scale_fill_gradient2(low="#5283ff",mid ="white", high ="#ff5c5c") +
      theme(axis.text.y = element_text(size  = 5)) +
      theme(legend.position = "bottom" )
    
  dev.off()
  
  
  # --------------- Check specific tissue marker --------------- #
   
  pdf(
    file = paste0(Save.Path,"/",ProjectName,"_CT/",ProjectName,"_nlDR_CTMarker.pdf"),
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
    file = paste0(Save.Path,"/",ProjectName,"_CT/",ProjectName,"_Heatmap_CellType_top",top_NSet,".pdf"),
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
    file = paste0(Save.Path,"/",ProjectName,"_CT/",ProjectName,"_Heatmap_CellType_top",top_NSet,".pdf"),
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
          pt.size =2) %>% BeautifyUMAP(FileName = paste0("/",Version,"/",ProjectName,"_CT/",ProjectName,"_nlDR_CellType"))
  DimPlot(scRNA.SeuObj,group.by = colnames(list_files.df)[2],  
          pt.size =0.5) %>% BeautifyUMAP(FileName = paste0("/",Version,"/",ProjectName,"_CT/",ProjectName,"_nlDR_",colnames(list_files.df)[2]))
  DimPlot(scRNA.SeuObj,group.by = "seurat_clusters",label.size = 7, label = TRUE,  
          pt.size =1) %>% BeautifyUMAP(FileName = paste0("/",Version,"/",ProjectName,"_Clusters/",ProjectName,"_nlDR_Clusters"))
  

  pdf(
    file = paste0(Save.Path,"/",ProjectName,"_CT/",ProjectName,"_nlDR_CellType_Sup.pdf"),
    width = 12,  height = 8
  )
  
    DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>% 
      BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)
    
    
    for (i in 1:(ncol(list_files.df)-1)) {
      print(DimPlot(scRNA.SeuObj, reduction = "umap", group.by = colnames(list_files.df)[i+1]) %>% 
              BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+ 
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2, split.by = colnames(list_files.df)[i+1], label = TRUE, label.size = 4) %>% 
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
      
      print(DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = colnames(list_files.df)[i+1]) %>% 
              BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.8, 0.15),AxisTitleSize=1.2, LegTextSize = 18)+ 
              theme(plot.title = element_text(vjust = 0.85)))
      print(DimPlot(scRNA.SeuObj, reduction = "tsne", ncol = 2, split.by = colnames(list_files.df)[i+1], label = TRUE, label.size = 4) %>% 
              BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                             SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9,OL_Thick = 1.5))
      
    }
    rm(i)
  
  dev.off()
  

  ## DotPlot
  DotPlot_Color1.set <- c("#de3767", "#de3767", "#4169e1", "#4169e1")
  DotPlot_Color2.set <- c("#5b8e7d","#7b2cbf")
  DotPlot_Color3.set <- c("#de3767", "#4169e1")
  
  pdf( 
    file = paste0(Save.Path,"/",ProjectName,"_CT/",ProjectName,"_DotPlot_CellType",".pdf"),
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
  ## Annotation Summary Table
    for (i in 1:(length(list_files.df)-1)) {
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
      Anno_Tar.set <- Anno.df[, ClassSet1] %>% unique()
      Anno_Tar.Num <- Anno_Tar.set %>% length()
      Anno_Tar_df.lt <- list()
      
      for (i in 1:Anno_Tar.Num) {
        Anno_Tar_df.lt[[i]] <- Anno.df[Anno.df[, ClassSet1]==Anno_Tar.set[i],]
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
          colnames(Anno_Freq_Tar_df.lt[[j]][[i]])[1] <- colnames(list_files.df)[ClassSet1+1]
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
      
    #### LinePlot ####
      # https://ithelp.ithome.com.tw/articles/10186047
      # Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
      #                                    levels = sort(unique(as.character(Freq_All.df$Cell_Type))))
  
        Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
                                        levels = Cell_Type_Order.set)
        
        ## Plot by State
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
      
      write.table( Freq_All_Cla.lt[["Freq_All_Cla.df"]] ,
                   file = paste0(Save.Path,"/",ProjectName,"_CellCount_CT_Cla.tsv"),
                   sep = "\t",
                   quote = F,
                   row.names = F
      )
      
      
      
      
      #### LinePlot ####
        Freq_All_Cla.df <- Freq_All_Cla.lt[["Freq_All_Cla.df"]]
        
        # Freq_All_Cla.df$Cell_Type <- factor(Freq_All_Cla.df$Cell_Type,
        #                                    levels = sort(unique(as.character(Freq_All_Cla.df$Cell_Type))))
        Freq_All_Cla.df$Cell_Type <- factor(Freq_All_Cla.df$Cell_Type,
                                           levels = Cell_Type_Order.set)
        
        CellNum_P3 <- ggplot(Freq_All_Cla.df, aes(x = factor(Cell_Type), y = Number, 
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
        
        CellNum_P4 <- ggplot(Freq_All_Cla.df, aes(x = factor(Cell_Type), y = Percent, 
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
      
        rm(Freq_All_Cla.df)
      
      
      ##### BarPlot ##### 
        # https://blog.gtwang.org/r/ggplot2-tutorial-layer-by-layer-plotting/3/
        colnames(Anno.df)[ncol(Anno.df)] <- "Cell_Type"
        Anno.df$Cell_Type <- factor(Anno.df$Cell_Type,
                                    levels = sort(unique(as.character(Anno.df$Cell_Type))))
        
        # sample
        BarPlot1_1 <- ggplot(Anno.df, aes(Cell_Type, fill=Anno.df[,1])) + 
          geom_bar(position="dodge")+theme_bw()+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
        BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                      XtextSize=18,  YtextSize=,18, xangle = 90,
                                      LegTextSize = 15)  + 
          theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
          labs(fill=colnames(Anno.df)[1]) -> BarPlot1_1
        BarPlot1_1
        
        ##### Export PDF file #####
        pdf(file = paste0(Save.Path,"/",ProjectName,"_CellCount_LinePlot.pdf"),
            width = 7, height = 7 )
          CellNum_P4
          CellNum_P3
          CellNum_P1
          CellNum_P2
          for (i in 1:(ncol(Anno.df)-1)) {
              # sample
              BarPlot1_1 <- ggplot(Anno.df, aes(Cell_Type, fill = Anno.df[,i])) + 
                geom_bar(position="dodge")+theme_bw()+
                theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
              BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                            XtextSize=18,  YtextSize=,18, xangle = 90,
                                            LegTextSize = 15)  + 
                theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
                labs(fill=colnames(Anno.df)[i]) -> BarPlot1_1
              print(BarPlot1_1)
              
              BarPlot1_2 <- ggplot(Anno.df, aes(Cell_Type, fill = Anno.df[,i])) + 
                geom_bar(position="fill")+theme_bw()+
                theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
              BarPlot1_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                            XtextSize=18,  YtextSize=,18, xangle = 90,
                                            LegTextSize = 15)  + 
                theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) +
                labs(fill=colnames(Anno.df)[i]) -> BarPlot1_2
              print(BarPlot1_2)
              rm(BarPlot1_1,BarPlot1_2)
              
            }
            rm(i)

        dev.off() # graphics.off()
        
        rm(CellNum_P1,CellNum_P2,CellNum_P3,CellNum_P4)
        
save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))

##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  #### Define group by different phenotype ####
    source("FUN_Find_Markers.R")
    scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
    
    # for (i in 1:(ncol(list_files.df)-1)) {
    #   scRNA.SeuObj[[paste0("celltype.",colnames(list_files.df)[1+i])]] <- paste(Idents(scRNA.SeuObj), as.matrix(scRNA.SeuObj[[colnames(list_files.df)[1+i]]]), sep = "_")
    #   
    # }
    
    scRNA.SeuObj[[paste0("celltype.",ClassSet2)]] <- paste(Idents(scRNA.SeuObj), 
                                                           as.matrix(scRNA.SeuObj[[ClassSet2]]), sep = "_")
    scRNA.SeuObj[[paste0("celltype.",ClassSet2,".",ClassSet3)]] <- paste(Idents(scRNA.SeuObj), 
                                                                         as.matrix(scRNA.SeuObj[[ClassSet2]]), 
                                                                         as.matrix(scRNA.SeuObj[[ClassSet3]]), sep = "_")
    Idents(scRNA.SeuObj) <- paste0("celltype.",ClassSet2,".",ClassSet3)

    
    DefaultAssay(scRNA.SeuObj) <- "RNA"
    
    
    ClassSet2.set <- list_files.df[[ClassSet2]] %>% unique()
    ClassSet3.set <- list_files.df[[ClassSet3]] %>% unique()
    CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
    # CellType.list <- CellType.list[-9]
    ####-------------- Find Marker gene in Male --------------####
      Sep_Cla3_FMar.Path <- paste0(Sampletype,"_",ProjectName,"_Separate_",ClassSet3.set[1],"_FindMarkers")
      dir.create(paste0(Save.Path,"/",Sep_Cla3_FMar.Path))
      
      # About 15 mins
      CCMarker_Male.lt <- list()
      for(i in c(1:length(CellType.list))){ 
        try({
          CCMarker_Male.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                                paste0(CellType.list[i],"_",ClassSet2.set[1],"_",ClassSet3.set[1]), 
                                                paste0(CellType.list[i],"_",ClassSet2.set[2],"_",ClassSet3.set[1]),
                                                CellType.list[i],
                                                Path = Save.Path,
                                                ResultFolder = paste0(Sep_Cla3_FMar.Path),
                                                ProjectTitle = ProjectName)
          # names(CCMarker_Male.lt)[[i]] <- paste0("CCMarker_Male.lt.",CellType.list[i])
          names(CCMarker_Male.lt)[[i]] <- paste0(CellType.list[i])
        })
      }
      rm(i)
      
      CCMarker_Male.lt <- CCMarker_Male.lt[!unlist(lapply(CCMarker_Male.lt,is.null))]
      
      
      ## Generate pdf and tif file for Male VolcanoPlot
      Sep_Cla3_Volcano.Path <- paste0(Sampletype,"_",ProjectName,"_Separate_",ClassSet3.set[1],"_VolcanoPlot")
      dir.create(paste0(Save.Path,"/",Sep_Cla3_Volcano.Path ))

      pdf(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),
          width = 7, height = 7 )
      for (i in 1:length(CellType.list)) {
        try({
          print(VolcanoPlot(CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S")]],
                            CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                            CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S_Neg_List")]], ShowGeneNum = 6)+ 
                  ggtitle(paste0(ProjectName,"_",ClassSet3.set[1],"_",CellType.list[i]))
          )
        })
      }
      dev.off() # graphics.off()
      rm(i)
      
      for (i in 1:length(CellType.list)) {
        try({
          tiff(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",CellType.list[i],".tif"), 
               width = 17, height = 17, units = "cm", res = 200)
          print(VolcanoPlot(CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S")]],
                            CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                            CCMarker_Male.lt[[i]][[paste0(ProjectName, "Marker.S_Neg_List")]])+ 
                  ggtitle(paste0(ProjectName,"_",ClassSet3.set[1],"_",CellType.list[i]))
          )
          
          graphics.off()
        })
      }
      rm(i,Sep_Cla3_FMar.Path)
    
    ####-------------- Find Marker gene in Female --------------####
      Sep_Cla3_FMar.Path <- paste0(Sampletype,"_",ProjectName,"_Separate_",ClassSet3.set[2],"_FindMarkers")
      dir.create(paste0(Save.Path, "/", Sep_Cla3_FMar.Path))
      
      # About 15 mins
      CCMarker_Female.lt <- list()
      
      for(i in c(1:length(CellType.list))){ 
        try({
          CCMarker_Female.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                                  paste0(CellType.list[i],"_",ClassSet2.set[1],"_",ClassSet3.set[2]), 
                                                  paste0(CellType.list[i],"_",ClassSet2.set[2],"_",ClassSet3.set[2]),
                                                  CellType.list[i],
                                                  Path = Save.Path,
                                                  ResultFolder = paste0(Sep_Cla3_FMar.Path),
                                                  ProjectTitle = ProjectName)
          # names(CCMarker_Female.lt)[[i]] <- paste0("CCMarker_Female.lt.",CellType.list[i])
          names(CCMarker_Female.lt)[[i]] <- paste0(CellType.list[i])
        })
      }
      rm(i)
      
      CCMarker_Female.lt <- CCMarker_Female.lt[!unlist(lapply(CCMarker_Female.lt,is.null))]
      
      
      ## Generate pdf and tif file for Female VolcanoPlot
      Sep_Cla3_Volcano.Path <- paste0(Sampletype,"_",ProjectName,"_Separate_",ClassSet3.set[2],"_VolcanoPlot")
      dir.create(paste0(Save.Path,"/",Sep_Cla3_Volcano.Path ))
      
      # dir.create(paste0(Save.Path,"/PBMC_SSA_Female_VolcanoPlot/"))
      
      pdf(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),
          width = 7, height = 7 )
        for (i in 1:length(CellType.list)) {
          try({
            print(VolcanoPlot(CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S")]],
                              CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                              CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S_Neg_List")]], ShowGeneNum = 6)+ 
                    ggtitle(paste0(ProjectName,"_",ClassSet3.set[2],"_",CellType.list[i]))
            )
          })
        }
        dev.off() # graphics.off()
        rm(i)
        
        for (i in 1:length(CellType.list)) {
          try({
            tiff(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",CellType.list[i],".tif"), 
                 width = 17, height = 17, units = "cm", res = 200)
            print(VolcanoPlot(CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S")]],
                              CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                              CCMarker_Female.lt[[i]][[paste0(ProjectName, "Marker.S_Neg_List")]])+ 
                    ggtitle(paste0(ProjectName,"_",ClassSet3.set[2],"_",CellType.list[i]))
            )
            
            graphics.off()
        })
      }
      rm(i,Sep_Cla3_FMar.Path)
      
save.image(paste0(Save.Path,"/08_1_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Separate).RData"))
    
    
##### 08_2 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ##-------------- Intersect_CellType --------------##
  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))
  
  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]
  
  CellType.list <- names(CCMarker_Male.lt)
  
  ##-------------- Venn Pos --------------##
  source("FUN_Venn.R")
  # pdf(file = paste0(Save.Path,"/PBMC_Female_VolcanoPlot.pdf"),width = 7, height = 7 )
  Sep_Cla3_Venn.Path <- paste0(Sampletype,"_",ProjectName,"_Separate_","_VennDiagrame")
  dir.create(paste0(Save.Path,"/",Sep_Cla3_Venn.Path))
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][[paste0(ProjectName, "Marker.S_Pos_List")]],
                                               CellType.list[i],"Pos","#9d0208","#f08080",SampleType = Sampletype,
                                               PathName = paste0(Save.Path,"/",Sep_Cla3_Venn.Path),
                                               ClassSet3_1 = ClassSet3.set[1], ClassSet3_2 =ClassSet3.set[2])
      names(Venn_CCMarker_Pos)[[i]] <- paste0("Venn_",ProjectName,"Marker.",CellType.list[i],"_Pos")
    })
  }
  rm(i)
  
  ##-------------- Venn Neg --------------##
  Venn_CCMarker_Neg <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      Venn_CCMarker_Neg[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][[paste0(ProjectName, "Marker.S_Neg_List")]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][[paste0(ProjectName, "Marker.S_Neg_List")]],
                                               CellType.list[i],"Neg","#00296b","#1368aa",SampleType=Sampletype,
                                               PathName = paste0(Save.Path,"/",Sep_Cla3_Venn.Path))
      
      names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_",ProjectName,"Marker.",CellType.list[i],"_Neg")
    })
  }
  rm(i,Sep_Cla3_Venn.Path)
  
  save.image(paste0(Save.Path,"/08_2_Find_",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_Venn.RData"))
  
  
##### 08_3 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  
  Idents(scRNA.SeuObj) <- "celltype.Cachexia"
  Sep_Cla3_FMar.Path <- paste0(Sampletype,"_",ProjectName,"_Pooled_FindMarkers")
  dir.create(paste0(Save.Path,"/",Sep_Cla3_FMar.Path))
  
  CCMarker_SPA.lt <- list()
  for(i in c(1:length(CellType.list))){ 
    try({
      CCMarker_SPA.lt[[i]] <- Find_Markers(scRNA.SeuObj, 
                                           paste0(CellType.list[i],"_",ClassSet2.set[1]), 
                                           paste0(CellType.list[i],"_",ClassSet2.set[2]),
                                           CellType.list[i],
                                           Path = Save.Path,
                                           ResultFolder =  Sep_Cla3_FMar.Path,
                                           ProjectTitle = ProjectName)
      
      # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
      names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i,Sep_Cla3_FMar.Path)
  
  CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]
  
  
  ## Generate pdf and tif file for VolcanoPlot
  Sep_Cla3_Volcano.Path <- paste0(Sampletype,"_",ProjectName,"_Pooled_","_VolcanoPlot")
  dir.create(paste0(Save.Path,"/",Sep_Cla3_Volcano.Path))
  pdf(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,".pdf"),width = 7, height = 7 )
    for (i in 1:length(CellType.list)) {
      try({
        print(VolcanoPlot(CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S")]],
                          CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Pos_List")]],
                          CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Neg_List")]], ShowGeneNum = 6)+ 
                ggtitle(paste0(Sampletype,"_",CellType.list[i]))
        )
      })
    }
    # graphics.off()
  dev.off()
  rm(i)
  
  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/",Sep_Cla3_Volcano.Path,"/",Sep_Cla3_Volcano.Path,"_",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Pos_List")]],
                        CCMarker_SPA.lt[[i]][[paste0(ProjectName,"Marker.S_Neg_List")]])+ 
                        ggtitle(paste0(Sampletype,"_",CellType.list[i]))
      )
      
      graphics.off()
    })
  }
  rm(i,Sep_Cla3_Volcano.Path)
  
save.image(paste0(Save.Path,"/08_3_Find__",Sampletype,"_",ProjectName,"marker_in_different_Cell_type_and_VolcanoPlot(Pooled).RData"))
  
################## (Pending) CCmarker matrix (Heatmap) ##################   
################## (Pending) CCmarker matrix LogFC (Heatmap) ##################  

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
  
  ##### Converting the Human gene name to Mouse gene name ##### 
    # #  Need to be optimized
    # # (Method1) bind the different length of column (Cannot use rbind)
    # # (Method2) Save the data as list first and than use do.call() to unlist to have dataframe
    # 
    # ## (Ori method)
    # Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*2))
    # for (i in 1:nrow(Pathway.all)) {
    #   #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    #   PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t() 
    #   colnames(PathwayN)="Temp"
    #   PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
    #   Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
    # }
    # 
    # Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
    # colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
    # 
    # rm(PathwayN)
    # 
    # # assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
    # # assign(colnames(Pathway.all)[i],Pathway.all[,i])
    
    ## (Method1)
    # Refer # https://stackoverflow.com/questions/3699405/how-to-cbind-or-rbind-different-lengths-vectors-without-repeating-the-elements-o
    # How to cbind or rbind different lengths vectors without repeating the elements of the shorter vectors?
    ## Modify by Charlene: Can use in MultRow
    bind_diff <- function(x, y){
      if(ncol(x) > ncol(y)){
        len_diff <- ncol(x) - ncol(y)
        y <- data.frame(y, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }else if(ncol(x) < ncol(y)){
        len_diff <- ncol(y) - ncol(x)
        x <- data.frame(x, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }
      rbind(x, y) 
    }
    
    ## Converting
    for (i in 1:nrow(Pathway.all)) {
      PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()  %>% as.data.frame()
      colnames(PathwayN)="Temp"
      PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
      PathwayN <- PathwayN[!PathwayN$MM.symbol  == 0,]
      PathwayN <- PathwayN[!is.na(PathwayN$MM.symbol),]
      if(i==1){
        Pathway.all.MM <- unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame()
      }else{
        Pathway.all.MM <- bind_diff(Pathway.all.MM,unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame())
        # Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
      }
    }
    
    Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
    colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
    #Pathway.all.MM[Pathway.all.MM==0] <-NA
    
    rm(PathwayN)
  
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
  
  
  
##### GO/Metascape #####  
  
  
  
##### inferCNV #####
  
##### Deconvolution #####
  
##### Beautify Figs #####  