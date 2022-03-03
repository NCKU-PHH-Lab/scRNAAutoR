# ## https://medium.com/datainpoint/r-essentials-writing-functions-f32d3c5cfc01
# # How to define a function
# FUNCTION_NAME <- function(INPUT1, INPUT2, ..., PARAM1, PARAM2, ...) {
#   # BODY
#   return(OUTPUT)
# }


#####----------------------------- Function -----------------------------#####

Venn_Intersect <- function(setA, setB,CellType.list,Type,ColorA,ColorB,SampleType="PBMC",PathName = getwd()) {
  ##-------------- Intersect --------------##
  Intersect_AB <- intersect(setA, setB)
  Union_AB <- union(setA, setB)
  Unique_A <-setdiff(setA, setB)
  Unique_B <-setdiff(setB,setA)
  
  Summary <- list(Intersect_AB,Union_AB,Unique_A,Unique_B)
  names(Summary) <- c('Intersect_AB','Union_AB','Unique_A','Unique_B')
  
  ##-------------- VennDiagram --------------##
  # Load the VennDiagram package
  # https://www.r-graph-gallery.com/14-venn-diagramm.html
  library(VennDiagram)
  
  # colorsT <- c("#2563e8", "#eb57c8")
  colorsT <- c(ColorA, ColorB)
  library("png")
  
  A <- na.omit(setA)
  B <- na.omit(setB)

   venn.diagram(x = list(A, B) ,
               category.names = c(paste0("M.",CellType.list), paste0("F.",CellType.list)),
               filename = paste0(PathName, "/IMGVenn_",SampleType,"_",Type,"_",CellType.list,".png"),
  #             output=TRUE,
               imagetype="png",
                height = 2400 , 
                width = 2400 , 
               scaled = FALSE,
               col = "Black",
                # Circles
                lwd = 2,
               fill = colorsT,

                # Numbers
                cex = 2,
                fontface = "bold",
                fontfamily = "sans",
  
               # Set names
           #    cat.col = 1,
               cat.col = colorsT,
               cat.cex = 1.5,
               cat.fontface = "bold",
               cat.fontfamily = "sans",

       #        cat.pos = c(90, 270),
               cat.default.pos = "outer",
               margin = 0.15, 
               cat.dist = c(0.15, 0.15),
               sub.just =c(1, 1)
  )

  
  # Display saved image
  options(repr.plot.height=12, repr.plot.width= 12)
  library("png")
  VennDiag <- readPNG(paste0(PathName,"/IMGVenn_",SampleType,"_",Type,"_",CellType.list,".png"))
  plot.new() 
  rasterImage(VennDiag,0,0,1,1)
  
  
  #####------------------------ Output ------------------------ #####
  OUTPUT <- list(Summary)
  names(OUTPUT) <- c('Summary')
  
  return(OUTPUT)
}

#####----------------------------- Function -----------------------------#####
