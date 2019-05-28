
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("msigdf"))
suppressPackageStartupMessages(library("openxlsx"))


# # conbind gene set have "_UP" "_DN" to "_ALL --------------------------------------------------------------------------------------------------------------
# hallmark <- msigdf.human %>% dplyr::select(geneset, symbol) %>% as.data.frame
# all_geneset <- hallmark$geneset %>% unique()
# up_geneset <- all_geneset[ endsWith(all_geneset, "_UP")] %>% substr(.,start=1,stop=nchar(.)-3 ) 
# dn_geneset <- all_geneset[ endsWith(all_geneset, "_DN")] %>% substr(.,start=1,stop=nchar(.)-3 ) 
# up_dn_geneset <- intersect(up_geneset,dn_geneset) 
# no_need_combine_geneset <- all_geneset[ all_geneset %in% up_dn_geneset ]
# need_combine_geneset <- up_dn_geneset[ ! up_dn_geneset %in% no_need_combine_geneset ]
# for (  rr in 1:length(need_combine_geneset)){
#   all_tbl <- hallmark %>% filter( geneset ==paste(need_combine_geneset[rr],"_UP", sep= "")| geneset ==paste(need_combine_geneset[rr],"_DN", sep= "") )  
#   all_tbl$geneset <- paste(need_combine_geneset[rr],"_ALL",sep = "" ) 
#   hallmark <- bind_rows(all_tbl, hallmark)
#   print(rr)
# } 
# -----------------------------------------------------------------------------------------------------------------------------------------------------------

option_list <- list( 
  make_option( c("-s","--samplePath")
              , type ="character"
              , help = "sample data(xlsx file)"),
  make_option( c("-g","--geneSetPath")
              , type ="character"
              , help = "gene set data(xlsx file)")
  )

load("/home/rstudio/gsea_msigdf_combine_up_down.Rdata")

opt <- parse_args(OptionParser(option_list=option_list))
  
if ( is.character(opt$samplePath) ){
  
  if(  is.character(opt$geneSetPath) ){
    gene_set <- read.xlsx(opt$geneSetPath)
    hallmark <- hallmark %>% as.data.frame %>% rbind(gene_set,.) 
  }
  

  samples <- read.xlsx(opt$samplePath) 

  DoGsea <- function(Values, GeneNames){
    names(Values) <- GeneNames
    Values <- Values[Values<1] 
    duplicate_check <- Values %>%.[duplicated(names(.))]
    
    if ( length(duplicate_check)!=0 ){
      return( paste("duplicared gene names", duplicate_check, sep = " ") )
    }
    Values <- sort(Values, decreasing = T)
    GSEA(Values, TERM2GENE = hallmark,pvalueCutoff = 0.1)
  }

  DoGseaKegg <- function(Values, GeneNames){
    names(Values) <- GeneNames
    Values <- Values[Values<1] 
    duplicate_check <- Values %>%.[duplicated(names(.))]
    
    if ( length(duplicate_check)!=0 ){
      return( paste("duplicared gene names", duplicate_check, sep = " ") )
    }
    
    toEntrezid <- data.frame(symbol = names(Values), values = Values) %>% 
      inner_join(.,  bitr( names(Values), fromType = "SYMBOL", toType = "ENTREZID",OrgDb="org.Hs.eg.db"),by= c("symbol"="SYMBOL") ) 
    
    Values <- toEntrezid$values 
    names(Values) <- toEntrezid$ENTREZID
    Values <- sort(Values, decreasing = T)
    
    gseKEGG(geneList     = Values,
            organism     = 'hsa',
            nPerm        = 1000,
            minGSSize    = 20,
            pvalueCutoff = 1,
            verbose      = FALSE) 
  }
 
  DoGseaGO <- function(Values, GeneNames, OntType){
    names(Values) <- GeneNames
    Values <- Values[Values<1] 
    duplicate_check <- Values %>%.[duplicated(names(.))]
    
    if ( length(duplicate_check)!=0 ){
      return( paste("duplicared gene names", duplicate_check, sep = " ") )
    }
    
    toEntrezid <- data.frame(symbol = names(Values), values = Values) %>% 
      inner_join(.,  bitr( names(Values), fromType = "SYMBOL", toType = "ENTREZID",OrgDb="org.Hs.eg.db"),by= c("symbol"="SYMBOL") ) 
    
    Values <- toEntrezid$values 
    names(Values) <- toEntrezid$ENTREZID
    Values <- sort(Values, decreasing = T)
    
    clusterProfiler::gseGO( geneList     = Values,
                            OrgDb        = org.Hs.eg.db,
                            ont          = OntType,
                            nPerm        = 1000,
                            minGSSize    = 100,
                            maxGSSize    = 500,
                            pvalueCutoff = 1,
                            verbose      = FALSE)
  }
  
  CombineRes <- function(Res){
    sample_name <- names(Res)
    fres <- Res[[1]]@result %>% mutate(sample = sample_name[1])
    for ( nn in 2:length(Res) ){
      fres <- rbind(fres, Res[[nn]]@result %>% mutate(sample = sample_name[nn]) ) 
    }
    fres %>% dplyr::select(sample,ID,Description,setSize,enrichmentScore,NES,pvalue,p.adjust,qvalues,rank,leading_edge,core_enrichment)
  }
  
 gsea_res <-  apply( samples[,-1],2,DoGsea, GeneNames = as.factor(samples[,1]) )
 kegg_res <-  apply( samples[,-1],2,DoGseaKegg, GeneNames = as.factor(samples[,1]) )
 go_cc_res <-  apply( samples[,-1],2,DoGseaGO, GeneNames = as.factor(samples[,1]), OntType = "CC" )
 go_bp_res <-  apply( samples[,-1],2,DoGseaGO, GeneNames = as.factor(samples[,1]), OntType = "BP" )
 go_mf_res <-  apply( samples[,-1],2,DoGseaGO, GeneNames = as.factor(samples[,1]), OntType = "MF" )
 
 gsea_fres <- CombineRes(Res = gsea_res)
 kegg_fres <- CombineRes(Res = kegg_res)
 go_cc_fres <- CombineRes(Res = go_cc_res)
 go_bp_fres <- CombineRes(Res = go_bp_res)
 go_mf_fres <- CombineRes(Res = go_mf_res)
 
 fres <- list( gsea_fres
               , kegg_fres
               , go_cc_fres
               , go_bp_fres
               , go_mf_fres)
 
 names(fres) <- c("GSEA","KEGG","GO Cellular component", "GO Biological process", "GO Molecular function")

 write.xlsx(fres, "gsea_result.xlsx")
 warnings()
} 
