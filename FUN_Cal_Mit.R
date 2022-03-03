scRNAMit <- function(PBMC.combined){

## https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md
## Calculating proportion of reads mapping to mitochondrial transcripts
library("AnnotationHub")

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)


# Check versions of databases available
ahDb %>% 
  mcols()


# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

#library("GPseq")
library(ensembldb)
# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")        


# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

#View(annotations)    


# # Extract IDs for mitochondrial genes
# mt <- annotations %>%
#   dplyr::filter(seq_name == "MT") %>%
#   dplyr::pull(gene_name)

# (Ch) Extract IDs for mitochondrial genes
mt.mm <- data.frame(gene = row.names(PBMC.combined@assays[["RNA"]]))
mt.mm <- mt.mm[grep("mt-",mt.mm$gene),]

# (Ch) Create metadata dataframe
metadata <- PBMC.combined@meta.data
# (Ch) Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# (Ch) Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

library(Matrix)
# # Number of UMIs assigned to mitochondrial genes
# metadata$mtUMI <- Matrix::colSums(counts[which(rownames(counts) %in% mt),], na.rm = T)

# (Ch) Number of UMIs assigned to mitochondrial genes
GeneExp.Mx <- as.matrix(PBMC.combined@assays[["RNA"]]@counts)
metadata$mtUMI <- Matrix::colSums(GeneExp.Mx[which(rownames(GeneExp.Mx) %in% mt.mm),], na.rm = T)

# Calculate of mitoRatio per cell
metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI


PBMC.combined@meta.data <- metadata

return(PBMC.combined)
}
