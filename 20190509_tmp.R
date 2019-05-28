
library(msigdf)
library(clusterProfiler)
library(tidyverse)
library(openxlsx)
library(optparse)



gene_set <- read.csv(file.choose(), header = F)

gene_set_test <-gene_set$V1 %>%as.character() %>% strsplit(.,split=", ",fixed=T) %>% unlist %>% gsub(" ", replacement = "_", .) %>% toupper() %>% gsub("/",replacement = "_",.) %>% gsub("-",replacement = "_",.)



test<-msigdf.human %>% filter(category_code == "hallmark")%>% select(geneset) %>% unique() %>%.$geneset%>% gsub("HALLMARK_",replacement = "",.)

have <- test[test %in% gene_set_test]
havnt<- gene_set_test[!gene_set_test %in% test] 
test[!test %in% gene_set_test]

sapply(havnt, function (x){test[grepl(x,test)]})


have <- c(have, "UV_RESPONSE_DN","UV_RESPONSE_UP","MYC_TARGETS_V1","MYC_TARGETS_V2","ESTROGEN_RESPONSE_LATE","ESTROGEN_RESPONSE_EARLY","TNFA_SIGNALING_VIA_NFKB","REACTIVE_OXIGEN_SPECIES_PATHWAY")
have = paste("HALLMARK_",have, sep = "")

data <- read.xlsx("22_24_29_34_P_vs_N_from_32pairs-without_OSCC4-Partek_YCY_OSCC_PDx_..._for_GSEA.xlsx",sheet = 1)


hallmark <- msigdf.human %>% filter(category_code == "hallmark" & geneset %in% have ) %>% select(geneset, symbol) %>% as.data.frame

data <- read.xlsx("22_24_29_34_P_vs_N_from_32pairs-without_OSCC4-Partek_YCY_OSCC_PDx_..._for_GSEA.xlsx",sheet = 1)
# remove duplicated gene name
#data <- dataOrigin %>% select(OSCC_22_P_vs_N, OSCC_24_P_vs_N, OSCC29_P_vs_N, OSCC34_P_vs_N,Gene.Symbol)
ptn22 <- data %>% filter(OSCC22_PT_P > 1 | OSCC22_PT_N >1) %>% mutate(genesymbol=Gene.Symbol, ratio=(OSCC22_PT_P/OSCC22_PT_N)) %>% select(genesymbol,ratio) %>% group_by(genesymbol) %>% filter(ratio==max(ratio))  %>% unique() %>% mutate(log2fc=log2(ratio)) %>% filter(log2fc!=0)
ptn24 <- data %>% filter(OSCC24_PT_P > 1 | OSCC24_PT_N >1) %>% mutate(genesymbol=Gene.Symbol, ratio=(OSCC24_PT_P/OSCC24_PT_N)) %>% select(genesymbol,ratio) %>% group_by(genesymbol) %>% filter(ratio==max(ratio))  %>% unique() %>% mutate(log2fc=log2(ratio)) %>% filter(log2fc!=0)
ptn29 <- data %>% filter(OSCC29_PT_P > 1 | OSCC29_PT_N >1) %>% mutate(genesymbol=Gene.Symbol, ratio=(OSCC29_PT_P/OSCC29_PT_N)) %>% select(genesymbol,ratio) %>% group_by(genesymbol) %>% filter(ratio==max(ratio))  %>% unique() %>% mutate(log2fc=log2(ratio)) %>% filter(log2fc!=0)
ptn34 <- data %>% filter(OSCC34_PT_P > 1 | OSCC34_PT_N >1) %>% mutate(genesymbol=Gene.Symbol, ratio=(OSCC34_PT_P/OSCC34_PT_N)) %>% select(genesymbol,ratio) %>% group_by(genesymbol) %>% filter(ratio==max(ratio))  %>% unique() %>% mutate(log2fc=log2(ratio)) %>% filter(log2fc!=0)

ptn22ByName <- ptn22$log2fc
ptn24ByName <- ptn24$log2fc
ptn29ByName <- ptn29$log2fc
ptn34ByName <- ptn34$log2fc

names(ptn22ByName) <- ptn22$genesymbol
names(ptn24ByName) <- ptn24$genesymbol
names(ptn29ByName) <- ptn29$genesymbol
names(ptn34ByName) <- ptn34$genesymbol

ptn22ByName <- sort(ptn22ByName, decreasing = TRUE)
ptn24ByName <- sort(ptn24ByName, decreasing = TRUE)
ptn29ByName <- sort(ptn29ByName, decreasing = TRUE)
ptn34ByName <- sort(ptn34ByName, decreasing = TRUE)

head(ptn22ByName)
head(ptn24ByName)
head(ptn29ByName)
head(ptn34ByName)

ptn22gseaMsigDbHallmark <- GSEA(ptn22ByName, TERM2GENE = hallmark,pvalueCutoff = 0.1)
ptn24gseaMsigDbHallmark <- GSEA(ptn24ByName, TERM2GENE = hallmark,pvalueCutoff = 0.1)
ptn29gseaMsigDbHallmark <- GSEA(ptn29ByName, TERM2GENE = hallmark,pvalueCutoff = 0.1)
ptn34gseaMsigDbHallmark <- GSEA(ptn34ByName, TERM2GENE = hallmark,pvalueCutoff = 0.1)

ptn22gseaMsigDbHallmark@result$sample="OSCC22"
ptn24gseaMsigDbHallmark@result$sample="OSCC24"
ptn29gseaMsigDbHallmark@result$sample="OSCC29"
ptn34gseaMsigDbHallmark@result$sample="OSCC34"

GSEA <- rbind(ptn22gseaMsigDbHallmark@result,
              ptn24gseaMsigDbHallmark@result,
              ptn29gseaMsigDbHallmark@result,
              ptn34gseaMsigDbHallmark@result)
GSEA 

GSEAPLOT <- GSEA %>% 
  dplyr::select(ID,NES,sample) %>%
  spread(sample,NES) %>%
  mutate(dummy1=NA) %>%
  gather("sample","NES",-ID)

ggplot(GSEAPLOT,aes(x=ID,y=sample,fill=NES))+
  geom_tile()+coord_polar()+
  scale_fill_gradient2(low="blue",mid="white",high="red",na.value = "white")+
  theme_bw()#+ theme(axis.text.y = element_blank(), axis.text.x = element_blank())



GSEAPLOT2 <- GSEA %>%   dplyr::select(ID,NES,sample) %>%  spread(sample,NES)


save(GSEA, file = "GSEA.Rdata")
