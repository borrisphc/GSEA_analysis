
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("openxlsx"))


option_list <- list(
  make_option(c("-f","--folderPath")
              , type = "character"
              , help = "Path of folder included sample file(csv file)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if ( is.character(opt$folderPath) ){
  
  setwd(as.character(opt$folderPath))
  
  files <-list.files()
  
  for( nn in 1:length(files)){
    setwd(as.character(opt$folderPath))
          gene_list <- read.csv(files[nn], header = F)
          
          gene_list <- gene_list$V1 %>% as.factor()
          
          
          used_ref_tbl <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID",OrgDb="org.Hs.eg.db")
          unused_ref_tbl <- data.frame(unused = gene_list[!gene_list %in% used_ref_tbl$SYMBOL])
          entrez_id <- used_ref_tbl$ENTREZID
          
          
          kegg_res <- enrichKEGG( gene = entrez_id
                                  , organism = "hsa"
                                  , pvalueCutoff = 0.01)
          
          all_entrez_id <- mappedkeys(org.Hs.egUNIGENE)
          
          BP <-
            clusterProfiler::enrichGO( gene          = entrez_id%>% as.numeric(),
                                       universe      = all_entrez_id,
                                       OrgDb         = 'org.Hs.eg.db',
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
          MF <-
            clusterProfiler::enrichGO( gene          = entrez_id%>% as.numeric(),
                                       universe      = all_entrez_id ,
                                       OrgDb         = 'org.Hs.eg.db',
                                       ont           = "MF",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
          
          CC <-
            clusterProfiler::enrichGO( gene          = entrez_id%>% as.numeric(),
                                       universe      = all_entrez_id ,
                                       OrgDb         = 'org.Hs.eg.db',
                                       ont           = "CC",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
          
          ALL <-
            clusterProfiler::enrichGO( gene          = entrez_id%>% as.numeric(),
                                       universe      = all_entrez_id ,
                                       OrgDb         = 'org.Hs.eg.db',
                                       ont           = "ALL",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
          
          
          if ( !is.null(kegg_res) ){
            for ( yy in 1:nrow(kegg_res@result) ){
              kegg_res@result$geneID[yy] <- kegg_res@result$geneID[yy] %>% strsplit(x = .,split = "/",fixed = T) %>% unlist %>% bitr(geneID = .,fromType = "ENTREZID", toType = "SYMBOL",OrgDb="org.Hs.eg.db") %>% .$SYMBOL %>% paste0(.,collapse = "/")
            }
          }
         
          
          check_null <- function (x){
            if (is.null(x)){c("no pathway enriched")}else{
              x@result
            }
          }
          
          fres <- list( used_ref_tbl
                        ,unused_ref_tbl
                        ,kegg_res %>% check_null
                        ,BP%>% check_null
                        ,MF%>% check_null
                        ,CC%>% check_null
                        ,ALL%>% check_null
                        )
         
          
          names(fres) <- c("used_gene","unused_gene","KEGG", "BP", "MF", "CC", "GO_ALL")
          setwd('/home/rstudio/')
          write.xlsx( fres, sub(pattern = ".csv",replacement = "_result.xlsx", files[nn]))
          print(nn)
          
  }
  
}







