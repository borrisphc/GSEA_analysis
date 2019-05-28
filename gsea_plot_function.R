# pptx. reference :https://www.brodrigues.co/blog/2018-10-05-ggplot2_purrr_officer/

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("officer"))
suppressPackageStartupMessages(library("rvg"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("pheatmap"))

option_list <- list( 
  make_option( c("-g","--GseaDataPath")
               , type ="character"
               , help = "Gsea result data(xlsx file)"))

opt <- parse_args(OptionParser(option_list=option_list))

if ( is.character(opt$GseaDataPath) ){
  
  Raw_data <- read.xlsx(opt$GseaDataPath) 
  
  create_pptx <- function(plot, path){
    if(!file.exists(path)) {
      out <- read_pptx()
    } else {
      out <- read_pptx(path)
    }
    
    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(plot), type = "body") %>% 
      print(target = path)
  }
  
  HeatmapDataParse <- function(DATA){
    HEATMAPPLOT <- DATA  %>% 
      dplyr::select(sample,ID,NES) %>%
      spread(sample,NES)
    rownames(HEATMAPPLOT) <- HEATMAPPLOT[,1]
    HEATMAPPLOT[is.na(HEATMAPPLOT)] <- 0
    
    HEATMAPPLOT[,-1]
  }

  heatmap_plot_data <- HeatmapDataParse(DATA = Raw_data)

  
  plot_res <- Raw_data  %>% 
    dplyr::select(sample,ID,NES) %>% 
    mutate(sample_plot = sample)%>%
    group_by(sample) %>% 
    nest %>% 
    mutate(plots = map2(.y = sample, .x = data, ~ ggplot(  data = .x %>%
                                                             spread(sample_plot,NES ) %>%
                                                             mutate(dummy1=NA) %>%
                                                             gather("sample_plot","NES",-ID) 
                                                           ,aes(x=ID,y=sample_plot,fill=NES) ) +
                          geom_tile()+coord_polar()+
                          scale_fill_gradient2(low="blue",mid="white",high="red",na.value = "white")+
                          theme_bw()
    ) #.map2    
    ) #mutate
  
  # read_pptx() %>%
  #   add_slide(layout = "Title and Content", master = "Office Theme") %>%
  #   ph_with_vg(code = print(heatmap_plot_res), type = "body") %>% 
  #   print(target ="/home/rstudio/heatmap_plot.pptx")
  pdf("heatmap.pdf", width=16+ncol(heatmap_plot_data), height=nrow(heatmap_plot_data)*0.05+13)
  pheatmap(heatmap_plot_data)
  graphics.off()
  map(plot_res$plots, create_pptx, path = "circle plot.pptx")
  
}

