options(repos = list(CRAN="http://cran.rstudio.com/"))
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
list.of.packages <- c(c("gridExtra","tidyverse"))
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}
library(tidyverse)
library(gridExtra)

#Genes to visualise
args <- commandArgs(trailingOnly = TRUE)
Gene_list <-  unlist(strsplit(args,","))
#Gene_list <-  c("NANOG","POU5F1")


#Function to plot genes
plot_markers <- function(df,markers,return_grid=T,type=c("tsne","pca","umap")){
  
  markers <- markers[markers %in% colnames(df)]
  plot_data_column = function (data, column){ ggplot(df)+geom_point(aes_string("tSNE1","tSNE2",colour= column))}
  myplots <- lapply(markers, plot_data_column, data = pca_df)
  # if T return a grid instead of list
  if(return_grid){
    n <- length(myplots)
    nCol <- floor(sqrt(n))
    grid<-do.call("grid.arrange", c(myplots, ncol=2))
    return(grid)
  }
  #return list 
  return(myplots)
}

Xiang_tsne <- read.csv("./Xiang_tsne.csv")
Xiang_tsne$lineage <- gsub("_Transition","",Xiang_tsne$lineage)
lineages <- c("eEPI","mEPI","lEPI","EPI-PrE","PrE","EPI-AMN","AMN","PGC","PSA","eTE","CTB","eSTB","STB","eEVT","EVT")
Xiang_tsne$lineage <- factor(Xiang_tsne$lineage,levels = lineages)
Xiang_tsne$date <- factor(Xiang_tsne$date, levels = paste0("D",sort(as.numeric(gsub("D","",unique(Xiang_tsne$date))))))

pdf("Boxplot.pdf",width = 15,height = 4*(length(Gene_list)/2))
Xiang_tsne %>% pivot_longer(A1BG:ZZZ3,names_to = "Gene",values_to = "FPKM") %>% 
  filter(Gene %in% Gene_list) %>% mutate(Gene = factor(Gene,levels = Gene_list)) %>% ggplot() +
  geom_boxplot(aes(x = lineage, y = FPKM, fill = lineage)) + facet_wrap(~Gene,scales = "free",ncol = 2,)
ggsave("Boxplot.pdf",width = 15,height = 4*(length(Gene_list)/2) ,dpi = 300)

pdf("tSNE.pdf",width = 15,height = 8*(length(Gene_list)/2))
plot_markers(Xiang_tsne,markers = c("lineage","date",Gene_list),type = "tsne")
dev.off()
