
# Libraries
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(topGO)

# Enrichment database
anno <- AnnotationDbi::select(org.Hs.eg.db,keys = keys(org.Hs.eg.db,keytype = "SYMBOL"),columns = c("ENSEMBL","SYMBOL","GENENAME"),keytype = "SYMBOL")
bg <- unique(anno$SYMBOL)

# Reactable functions
render.reactable.cell.with.tippy <- function(text,tooltip){
  
  div(
    style = "text-decoration: underline;
                text-decoration-style: dotted;
                text-decoration-color: #FF6B00;
                cursor: info;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;",
    tippy(text = text,tooltip = tooltip)
  )
  
}

dataListFilter <- function(tableId,style = "width: 100%; height: 28px;") {
  function(values,name) {
    
    dataListId <- sprintf("%s-%s-list",tableId,name)
    tagList(
      tags$input(
        type = "text",
        list = dataListId,
        oninput = sprintf("Reactable.setFilter('%s', '%s', event.target.value || undefined)",tableId,name),
        "aria-label" = sprintf("Filter %s",name),
        style = style
      ),
      tags$datalist(
        id = dataListId,
        lapply(unique(values),function(value) tags$option(value = value))
      )
    )
  }
  
}

bar_chart <- function(label,width = "100%",height = "1rem",fill = "#00bfc4",background = NULL) {
  
  bar <- div(
    style = list(background = fill,width = width,height = height)
  )
  chart <- div(
    style = list(flexGrow = 1,marginLeft = "0.5rem",background = background),
    bar
  )
  div(
    style = list(display = "flex",alignItems = "center"),
    label,
    chart
  )
  
}

# GO enrichment function
Run_enrichment <- function(sig_genes,anno) {
  
  an_sig <- as.data.frame(subset(anno,SYMBOL %in% sig_genes))
  in_universe <- bg %in% c(an_sig$SYMBOL,bg) 
  in_selection <-  bg %in% an_sig$SYMBOL 
  alg <- factor(as.integer(in_selection[in_universe]))
  names(alg) <- bg[in_universe]
  tgd <- new("topGOdata",ontology = "BP",allGenes = alg,nodeSize = 5,annot = annFUN.org,mapping = "org.Mm.eg.db",ID = "symbol")
  result_topgo <- runTest(tgd,algorithm = "classic",statistic = "Fisher")
  all_go <- usedGO(tgd)
  go_results <- GenTable(tgd,Fisher.classic = result_topgo,orderBy = "Fisher.classic",topNodes = length(all_go),numChar = 10000)
  go_results$Fisher.classic <- ifelse(grepl("<",go_results$Fisher.classic),0,go_results$Fisher.classic)
  go_results$Fisher.classic <- round(p.adjust(go_results$Fisher.classic,method = "BH"),digits = 4)
  go_results <- go_results[go_results$Fisher.classic < 0.1,]
  names(go_results) <- c("Term Id","Term label","# Reference","# Genes","Fold enrichment","Adj. p-value")
  go_results$`Adj. p-value` <- round(go_results$`Adj. p-value`,3)
  go_results$`Fold enrichment` <- round(go_results$`Fold enrichment`,3)
  return(go_results)
  
}

