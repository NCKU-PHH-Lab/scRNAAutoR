BeautifyUMAP = function( P1,FileName = "NO1",RVer = Version,
                           LegPos = c(1.01, 0.5),LegBox = "vertical",LegDir="vertical",
                           LegTextSize = 22,
                           TH= 0.5,TV= 0.5, TitleSize = 30,
                           XtextSize=25,  YtextSize=25,  XaThick=0.9,  YaThick=0.9, xangle =0,
                           AxisTitleSize=1.5, AspRat=1,SubTitSize = 15 
                           ){
  
  library(ggplot2)
  library(graphics)
  
  P2 <-  P1 +            
    theme(axis.text.x = element_text(size = XtextSize,angle = xangle, hjust = 1, vjust = .5), #face="bold", # Change the size along the x axis
          axis.text.y = element_text(size = YtextSize), # Change the size along the y axis
          
          axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"), 
          axis.title = element_text(size = rel(AxisTitleSize),face="bold"),
          plot.title = element_text(color="black", 
                                    size=TitleSize, 
                                    face="bold.italic",
                                    hjust = TH,vjust =TV), # margin = margin(t = 0.5, b = -7),
          #     plot.background = element_rect(fill = 'chartreuse'),
          legend.title = element_text(size=12, color = "black", face="bold"),
          legend.text = element_text(colour="black", size= LegTextSize,face="bold"),
          legend.background = element_rect(fill = alpha("white", 0.5)),
          #      legend.position = c(0.1, 0.18),
          #     plot.text = element_text(size = 20),
          aspect.ratio=AspRat) + #square plot
    theme(axis.line.x = element_line(colour = "black", size = XaThick),
          axis.line.y = element_line(colour = "black", size = YaThick))+
    theme(legend.position = LegPos , legend.box = LegBox ,legend.direction=LegDir)+ # vertical, horizontal
    #theme(legend.position = "top" , legend.direction = LegDir) # vertical, horizontal
    theme(strip.text = element_text(size=SubTitSize))+ 
    theme(panel.background = element_rect(colour = "black",size = 1.5))

  ##
  P3 <- P1 + theme(panel.grid.major =element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),#去除背景
                    panel.border = element_blank(), #去除邊框
                   aspect.ratio=AspRat)+
              theme(axis.line.x = element_line(colour = "black", size = 0.8),
                    axis.line.y = element_line(colour = "black", size = 0.8),
                    panel.background = element_rect(colour = "black",size = 1.5),
                    axis.text.x = element_text(size = 0),
                    axis.text.y = element_text(size = 0),
                    axis.title = element_text(size = 0),
                    plot.title = element_text(size=0))+
              theme(legend.position="none")+
              theme(axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())
    
  
  ## Export PDF
  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/UMAP",FileName,".pdf"),
    width = 12,  height = 8
  )
  print(P2)
  print(P3)
  dev.off()
  # graphics.off()
  
  ## Export tiff
  tiff(
    file = paste0(setwd(getwd()),"/",Version,"/UMAP",FileName,"01.tiff"), 
    width = 28, height = 17, units = "cm", res = 200)
    print(P2)
   graphics.off()
  tiff(
    file = paste0(setwd(getwd()),"/",Version,"/UMAP",FileName,"02.tiff"), 
    width = 17, height = 17, units = "cm", res = 200)
   print(P3)
   graphics.off()
  
  return(P2)
}
