
# To do
# fix edge weight colors
# summary of number of up/dpwn regulated?
# readme for different types of datasets (how to best utilize the app for your usecase)
  # possibly use a vignette
# start a github repo for this project

# library(BiocManager)
# options(repos = BiocManager::repositories())

# Libraries
library(shiny)
library(shinyWidgets)
library(RColorBrewer)
library(sparkline)
library(bslib)
library(scales)
library(ggplot2)
library(plotly)
library(dplyr)
library(networkD3)
library(httr)
library(jsonlite)
library(webshot2)
library(htmlwidgets)
library(htmltools)
library(tippy)
library(reactable)
library(BiocManager)

# Options
options(shiny.maxRequestSize = 100*1024^2)

# Load in scripts
source("utils.R") # /Users/8yt/Documents/Network Application/

# Requirements
# (1) Must make sure that the input file has column names "Group" and "Gene" corresponding to omics type and gene name
# (2) Must make sure that the input file is a tab separated .txt file
# (4) For the multi-omics tab data make sure the input file has unique rows only

# Counter for tabs
tab_count <- 0
prev_tabs <- list()


# Define UI 
ui <- fluidPage(
  
    theme = bs_theme(version = 4,bootswatch = "minty",primary = "#004EFF", secondary = "#004EFF"),
    
    chooseSliderSkin("Flat",color = "#004EFF"),
    
    navbarPage("Network Application",
      
      tabPanel("Multi-omic Integration",
               
      tags$head(tags$style(".btn { width: 100%; background-color: #004EFF; }")),
      tags$head(tags$style(HTML("div#inline label { width: 25%; font-size: 12px; }"))),
      
      sidebarLayout(
          
          sidebarPanel(
              fileInput(
                  "file1",
                  "Upload file",
                  accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv"
                  )
              ),
              uiOutput("factor_filters"),
              tabsetPanel(
                  id = "data_filters",
                  type = "pills",
                  selected = NULL
              ),
              uiOutput("gene_filter"),
              fluidRow(
                  column(6,downloadButton('download_network','network')),
                  column(6,downloadButton('download_data','data'))
              ),
              width = 4
          ),
  
          mainPanel(
              forceNetworkOutput("net_visual",width = "100%",height = "650px"),
              reactableOutput("network_table",width = "100%"),
              width = 8
          )
        
        )
      
      ),
      
      tabPanel("Gene Enrichment",
               
      tags$head(tags$style(".btn { width: 100%; }")),
      tags$head(tags$style(HTML("div#inline label { width: 25%; font-size: 12px; }"))),
      
      sidebarLayout(
         
          sidebarPanel(
             uiOutput("factor_filters_go"),
             tabsetPanel(
               id = "data_filters_go",
               type = "pills",
               selected = NULL
             ),
             actionButton("go","run enrichment"),
             br(),br(),
             downloadButton('download_enrichment','enrichment results'),
             width = 4
          ),
        
          mainPanel(
            reactableOutput("go_table",width = "100%"),
            width = 8
          )
      
      )
    
    )
  
  )
  
)

# Define server logic required to draw a histogram
server <- function(input,output,session) {
  
    data <- reactive({
        
        inFile <- input$file1
        if (is.null(inFile))
            return(NULL)
        if(grepl(".txt",inFile$datapath)) {
          sep <- "\t"
        } else {
          sep <- ","
        }
        data <- unique(read.table(inFile$datapath,header = TRUE,sep = sep,stringsAsFactors = TRUE))
        if(sum(grepl("ENSG\\d+",data$Gene)) == 0) {
            data$Gene.Symbol <- data$Gene
        } else {
            r <- POST(url,body = list(api = 1,ids = toJSON(data$Gene)))
            output <- fromJSON(httr::content(r,"text"),flatten = TRUE)
            Symbols <- data.frame("ID" = names(output),"Gene Symbol" = Reduce(c,output))
            data <- merge(
                data,
                Symbols,
                by.x = "Gene",
                by.y = "ID",
                all.x = TRUE
            )
        }
        data$Gene.Symbol <- ifelse(is.na(as.character(data$Gene.Symbol)),as.character(data$Gene),as.character(data$Gene.Symbol))
        data <- data[,c("Group","Gene","Gene.Symbol",names(data)[!(names(data) %in% c("Group","Gene","Gene.Symbol"))])]
        groups <- unique(sapply(strsplit(unique(as.character(data$Group)),"_"),"[[",1))
        group_names <- case_when(
            grepl("RNA",groups,ignore.case = TRUE) ~ "RNA",
            grepl("CHIP",groups,ignore.case = TRUE) ~ "ChIP",
            grepl("ATAC",groups,ignore.case = TRUE) ~ "ATAC",
            grepl("DNA-meth",groups,ignore.case = TRUE) ~ "DNA",
            grepl("GWAS",groups,ignore.case = TRUE) ~ "GWAS",
            grepl("LC",groups,ignore.case = TRUE) ~ "PRO"
        )
        columns <- names(data)[!(names(data) %in% c("Group","Gene","Gene.Symbol"))]
        factors <- columns[do.call("c",lapply(columns,function(x){is.factor(data[[x]])}))]
        data$Label <- paste0(data$Group," (",apply(data.frame(data[,factors]),1,function(x){
          paste0(factors,": ",x,collapse = ", ")
        }),")")
        data$Label_ <- apply(data[,c("Group",factors)],1,paste0,collapse = "_")
        values <- columns[do.call("c",lapply(columns,function(x){is.numeric(data[[x]])}))]
        values_select <- lapply(1:length(groups),function(group){
            values_select <- data.frame(data[grepl(groups[group],data$Group),values])
            values[apply(values_select,2,function(col){all(!is.na(col))})]
        })
        value_names <- lapply(1:length(values_select),function(value){
            case_when(
                grepl("log2FC",values_select[[value]],ignore.case = TRUE) ~ "Log<sub>2</sub>(FC)&emsp;"
                ,
                grepl("adj.pval|adjpval|adjp|pvaladj|padj",values_select[[value]],ignore.case = TRUE) ~ "Adj. p-value&emsp;"
            )
        })
        # fix underlying dataset later
        # if("Cell.Type" %in% names(data)) {
        #   data$Cell.Type <- gsub(" [0-9]+","",data$Cell.Type)
        #   data$Cell.Type <- as.factor(data$Cell.Type)  
        # }
        data <- list(data,groups,values_select,group_names,value_names,factors) # ,tab_count,prev_tabs
        names(data) <- c("data","groups","values_select","group_names","value_names","factors") # ,"tab_count","prev_tabs"
        return(data)
        
    })
    
    output$factor_filters <- renderUI({
        
        req(!is.null(input$file1))
        data <- data()$data
        factors <- c("Group",data()$factors)
        lapply(factors,function(factor){
            selectInput(
                inputId = factor,
                label = gsub("_"," ",factor),
                choices = as.character(unique(data[[factor]])),
                selected = unique(data[[factor]])[1],
                multiple = TRUE
            )    
        })
        
    })
    
    observeEvent(input$file1,{

        req(!is.null(input$file1))
        prev_tabs[[as.character(tab_count)]] <<- data()$group_names
        tab_count <<- tab_count + 1
        data <- data()$data
        groups <- data()$groups
        group_names <- data()$group_names
        values_select <- data()$values_select
        value_names <- data()$value_names
        if(tab_count > 1) {
          lapply(1:length(prev_tabs[[as.character(tab_count - 2)]]),function(tab){
            removeTab(
              inputId = "data_filters",
              target = prev_tabs[[as.character(tab_count - 2)]][tab]
            )
          })
        }
        lapply(1:length(groups),function(group){
          appendTab(
            inputId = "data_filters",
            tabPanel(
              group_names[group],
              lapply(1:length(value_names[[group]]),function(value){
                if(grepl("Adj. p-value",value_names[[group]][value])) {
                  div(id = "inline",
                      selectInput(
                        inputId = paste0(groups[group],"_",values_select[[group]][value]),
                        label = HTML(value_names[[group]][value]),
                        choices = c(0.001,0.005,0.01,0.05,0.1),
                        selected = 0.05,
                        width = "100%"
                      )
                  )
                } else {
                  div(id = "inline",
                      sliderInput(
                        inputId = paste0(groups[group],"_",values_select[[group]][value]),
                        label = HTML(value_names[[group]][value]),
                        min = round(min(data[grepl(groups[group],data$Group),values_select[[group]][value]]),2),
                        max = round(max(data[grepl(groups[group],data$Group),values_select[[group]][value]]),2),
                        value = c(
                          -0.1,
                          0.1
                        ),
                        width = "100%"
                      )
                  )
                }
              })
            ),select = TRUE
          )
        })

    })

    output$gene_filter <- renderUI({

      req(!is.null(input$file1))
      textInput(
        inputId = "gene_filter",
        label = "Explore gene(s):",
        placeholder = ""
      )

    })
    
    data_filtered <- reactive({
        
        req(!is.null(input$file1))
        data <- data()$data
        groups <- data()$groups
        values_select <- data()$values_select
        value_names <- data()$value_names
        factors <- c("Group",data()$factors)
        labels <- apply(expand.grid(lapply(factors,function(x){input[[x]]})),1,paste0,collapse = "_")
        data <- data %>% filter(Label_ %in% labels) %>% data.frame()
        data_filtered <- list()
        for(group in 1:length(groups)) {
          data_subset <- data[data$Group == groups[group],]
            for(value in 1:length(value_names[[group]])) {
                if(!grepl("Adj. p-value",value_names[[group]][value])) {
                  data_subset <- data_subset[
                      data_subset[[values_select[[group]][value]]] <= input[[paste0(groups[group],"_",values_select[[group]][value])]][1] |
                      data_subset[[values_select[[group]][value]]] >= input[[paste0(groups[group],"_",values_select[[group]][value])]][2]
                    ,]
                } else {
                  data_subset <- data_subset[
                      data_subset[[values_select[[group]][value]]] <= as.numeric(input[[paste0(groups[group],"_",values_select[[group]][value])]])
                    ,]
                }
            }
            data_filtered[[group]] <- data_subset
        }
        data_filtered <- do.call("rbind",data_filtered)
        if(any(grepl("log2fc",names(data_filtered),ignore.case = TRUE))) {
          data_filtered$Value <- data_filtered[[names(data_filtered)[which(grepl("log2fc",names(data_filtered),ignore.case = TRUE))]]]
        } else {
          data_filtered$Value <- data_filtered[[names(data_filtered)[which(grepl("padj|adj.pval|adjp",names(data_filtered),ignore.case = TRUE))]]]
        }
        data_filtered <- data_filtered[,c("Label",factors,"Gene.Symbol","Value")]
        data_filtered
        
    })
    
    data_net <- reactive({
        
        req(input$file1)
        shiny::validate(need(nrow(data_filtered()) > 0,"No data for network. Adjust filters."))
        factors <- data()$factors
        data_filtered <- data_filtered()
        nodes <- data_filtered[,c("Label","Gene.Symbol","Value")]
        nodes$Type <- "gene"
        nodes <- do.call("rbind",lapply(unique(nodes$Gene.Symbol),function(x){
            g_ <- nodes[nodes$Gene.Symbol == x,]
            if(nrow(g_) == 1) {
                node_name <- paste0(g_$Gene.Symbol,": ",round(g_$Value,3))
                node <- g_$Gene.Symbol
                gs <- g_$Gene.Symbol
            } else {
                node_name <- paste0(unique(g_$Gene.Symbol),": ",paste0(paste0(g_$Label," [",round(g_$Value,3),"]"),collapse = ", "))
                node <- unique(g_$Gene.Symbol)
                gs <- unique(g_$Gene.Symbol)
            }
            data.frame("Node" = node,"Gene Symbol" = gs,"Node_name" = node_name,"Type" = "Gene")
        }))
        extra_nodes <- data.frame("Node" = unique(data_filtered$Label),"Gene.Symbol" = unique(data_filtered$Label),"Node_name" = unique(data_filtered$Label),"Type" = unique(data_filtered$Label))
        group <- do.call("c",lapply(strsplit(extra_nodes$Type," \\("),function(x){x[1]}))
        types <- strsplit(gsub(")","",gsub("\\w+-\\w+ \\(","",extra_nodes$Type)),", ")
        types <- lapply(types,function(type){gsub("\\w+[\\.]?\\w+: ","",type)})
        types <- do.call("rbind",lapply(1:length(types),function(x){
          dat <- data.frame(matrix(types[[x]],1))
          dat
        }))
        types <- data.frame(types[,do.call("c",lapply(1:ncol(types),function(x){length(unique(types[,x])) != 1}))])
        extra_nodes$Type <- paste0(group," (",apply(types,1,paste0,collapse = ", "),")")
        nodes <- rbind(nodes,extra_nodes)
        if(input$gene_filter != "") {
          explore <- do.call("c",strsplit(input$gene_filter," "))
          explore <- trimws(explore,"both")
          explore <- explore[which(explore != "")]
          explore <- paste0(explore,collapse = "|")
          if(any(grepl(explore,nodes$Gene.Symbol,ignore.case = TRUE))) {
            nodes[which(grepl(explore,nodes$Gene.Symbol,ignore.case = TRUE)),]$Type <- "Explore"
            colors <- brewer.pal(n = length(unique(nodes$Type)),name = 'Set2')
            colors[which(unique(nodes$Type) == "Explore")] <- "#FF0000"
            colors[which(unique(nodes$Type) == "Gene")] <- "#2374BA"
            color_string <- paste0("d3.scaleOrdinal().domain([",paste0("'",paste0(unique(nodes$Type),sep = "",collapse = "','"),"'"),"]).range([",paste0("'",paste0(colors,sep = "",collapse = "','"),"'"),"])")
          } else {
            colors <- brewer.pal(n = length(unique(nodes$Type)),name = 'Set2')
            colors[which(unique(nodes$Type) == "Gene")] <- "#2374BA"
            color_string <- paste0("d3.scaleOrdinal().domain([",paste0("'",paste0(unique(nodes$Type),sep = "",collapse = "','"),"'"),"]).range([",paste0("'",paste0(colors,sep = "",collapse = "','"),"'"),"])")
          }
        } else {
          colors <- brewer.pal(n = length(unique(nodes$Type)),name = 'Set2')
          colors[which(unique(nodes$Type) == "Gene")] <- "#2374BA"
          color_string <- paste0("d3.scaleOrdinal().domain([",paste0("'",paste0(unique(nodes$Type),sep = "",collapse = "','"),"'"),"]).range([",paste0("'",paste0(colors,sep = "",collapse = "','"),"'"),"])")
        }
        nodes$size <- ifelse(grepl("Gene|Explore",nodes$Type),3,70)
        nodes$`Node ID` <- 0:(nrow(nodes)-1)
        edges <- data_filtered[,c("Label","Gene.Symbol","Value")]
        names(edges) <- c("Source","Target","Value")
        edge_names <- names(edges)
        edges <- merge(
            edges,
            nodes,
            by.x = "Source",
            by.y = "Node",
            all.x = TRUE
        ) %>% dplyr::select(all_of(edge_names),"Node ID") %>% data.frame()
        names(edges)[ncol(edges)] <- "Source ID"
        edge_names <- names(edges)
        edges <- merge(
            edges,
            nodes[,c("Node","Node_name","Type","size","Node ID")],
            by.x = "Target",
            by.y = "Node",
            all.x = TRUE
        ) %>% dplyr::select(all_of(edge_names),"Node ID")
        names(edges)[ncol(edges)] <- "Target ID"
        edges$Link <- ifelse(grepl("RNA-seq",edges$Source) & all(edges$Value != 1),rescale(abs(edges$Value),to = c(1,30)),1)
        return <- list(data_filtered,nodes,edges,color_string)
        names(return) <- c("data_filtered","nodes","edges","color_string")
        return
                
    })
    
    download_network <- reactive({
        
        req(data_net())
        data <- data_net()
        nodes <- data$nodes
        edges <- data$edges
        link_color <- rescale(edges$Value)
        trans <- div_gradient_pal(low = "red",mid = "grey",high = "blue",space = "Lab")
        link_color <- trans(link_color)
        color_string <- data$color_string
        fn <- forceNetwork(
            Links = edges,
            Nodes = nodes,
            fontFamily = "Arial",
            Nodesize = "size",
            Source = "Source ID",
            Target = "Target ID",
            Value = "Link",
            linkColour = link_color,
            NodeID = "Node_name",
            charge = -100,
            Group = "Type",
            zoom = TRUE,
            opacity = 1,
            legend = TRUE,
            opacityNoHover = 0,
            colourScale = networkD3::JS(color_string)
        )
        fn <- htmlwidgets::onRender(fn,jsCode = '
            function (el, x) {
                d3.select("svg").append("g").attr("id","legend-layer");
                var legend_layer = d3.select("#legend-layer");
                d3.selectAll(".legend")
                .each(function() { legend_layer.append(() => this); });
            }
        ')
        fn
        
    })
    
    output$net_visual <- renderForceNetwork({
        
        req(download_network())
        download_network()
        
    })
    
    output$download_network <- downloadHandler(
        filename = function() {
            'network.html'
            # 'network.svg'
        },
        content = function(file) {
            saveNetwork(download_network(),file)
        }
    )
    
    output$network_table <- renderReactable({
        
        req(data_filtered())
        data <- data_filtered()
        data <- data[,names(data)[!(names(data) %in% c("Label","Label_"))]]
        data$Value <- round(data$Value,2)
        reactable(
            data,
            highlight = TRUE,
            bordered = TRUE,
            compact = TRUE,
            filterable = TRUE,
            showSortable = TRUE,
            defaultColDef = colDef(
              footer = function(values_select) {
                if (!is.numeric(values_select)) return()
                sparkline(values_select,type = "box",width = 100,height = 30)
              }
            ),
            columns = list(
              Value = colDef(
                name = "Log2FC",
                filterable = TRUE,
                filterMethod = JS(
                  "function(rows, columnId, filterValue) {
                    return rows.filter(function(row) {
                      return row.values[columnId] >= filterValue
                    })
                  }"
                )
              )
            )
        )
        
    })
    
    output$download_data <- downloadHandler(
        
        filename = function() {
            'network_data.txt'
        },
        content = function(file) {
            write.table(data_filtered(),file,sep = "\t",row.names = FALSE,col.names = TRUE)
        }
    
    )
    
    output$factor_filters_go <- renderUI({
      
      req(!is.null(input$file1))
      data <- data()$data
      factors <- c("Group",data()$factors)
      lapply(factors,function(factor){
        selectInput(
          inputId = paste0(factor,"_go"),
          label = gsub("_"," ",factor),
          choices = as.character(unique(data[[factor]])),
          selected = unique(data[[factor]])[1],
          multiple = FALSE
        )    
      })
      
    })
    
    observeEvent(input$file1,{
      
      req(!is.null(input$file1))
      data <- data()$data
      groups <- data()$groups
      group_names <- data()$group_names
      values_select <- data()$values_select
      value_names <- data()$value_names
      if(tab_count > 1) {
        lapply(1:length(prev_tabs[[as.character(tab_count - 2)]]),function(tab){
          removeTab(
            inputId = "data_filters_go",
            target = prev_tabs[[as.character(tab_count - 2)]][tab]
          )
        })
      }
      lapply(1:length(groups),function(group){
        appendTab(
          inputId = "data_filters_go",
          tabPanel(
            group_names[group],
            lapply(1:length(value_names[[group]]),function(value){
              if(grepl("Adj. p-value",value_names[[group]][value])) {
                div(id = "inline",
                    selectInput(
                      inputId = paste0(groups[group],"_",values_select[[group]][value],"_go"),
                      label = HTML(value_names[[group]][value]),
                      choices = c(0.001,0.005,0.01,0.05,0.1),
                      selected = 0.05,
                      width = "100%"
                    )
                )
              } else {
                div(id = "inline",
                    sliderInput(
                      inputId = paste0(groups[group],"_",values_select[[group]][value],"_go"),
                      label = HTML(value_names[[group]][value]),
                      min = round(min(data[grepl(groups[group],data$Group),values_select[[group]][value]]),2),
                      max = round(max(data[grepl(groups[group],data$Group),values_select[[group]][value]]),2),
                      value = c(
                        -0.1,
                        0.1
                      ),
                      width = "100%"
                    )
                )
              }
            })
          ),select = TRUE
        )
      })
      
    })
    
    go_results <- eventReactive(input$go,{
      
      # showModal(modalDialog("Running GO enrichment analysis. Please be patient.",footer = NULL,easyClose = TRUE,size = "m"))
      req(!is.null(input$file1))
      data <- data()$data
      groups <- data()$groups
      group_names <- data()$group_names
      values_select <- data()$values_select
      value_names <- data()$value_names
      factors <- c("Group",data()$factors)
      data_filtered <- list()
      for(factor in 1:length(factors)) {
        if(factor == 1) {
          data_filtered <- data[data[[factors[factor]]] == input[[paste0(factors[factor],"_go")]],]  
        } else {
          data_filtered <- data_filtered[data_filtered[[factors[factor]]] == input[[paste0(factors[factor],"_go")]],]
        }
      }
      for(value in 1:length(values_select[[which(groups == input$Group_go)]])) {
        if(!grepl("pval|padj|adjpval|adj.pval",values_select[[which(groups == input$Group_go)]][value],ignore.case = TRUE)) {
          data_filtered <- data_filtered[
            data_filtered[[values_select[[which(groups == input$Group_go)]][value]]] <= input[[paste0(input$Group_go,"_",values_select[[which(groups == input$Group_go)]][value],"_go")]][1] |
              data_filtered[[values_select[[which(groups == input$Group_go)]][value]]] >= input[[paste0(input$Group_go,"_",values_select[[which(groups == input$Group_go)]][value],"_go")]][2]
            ,]  
        } else {
          data_filtered <- data_filtered[
            data_filtered[[values_select[[which(groups == input$Group_go)]][value]]] <= input[[paste0(input$Group_go,"_",values_select[[which(groups == input$Group_go)]][value],"_go")]]
            ,] 
        }
      }
      shiny::validate(need(nrow(data_filtered) > 0,"No genes included for testing. Try adjusting the numeric thresholds."))
      sig_genes <- data_filtered$Gene.Symbol
      go_results <- Run_enrichment(sig_genes,anno)
      shiny::validate(need(nrow(go_results) > 0,"No significant GO terms. Try adjusting the numeric thresholds."))
      go_results$Genes <- do.call("c",lapply(1:nrow(go_results),function(x){
        paste0(intersect(sig_genes,genesInTerm(object = tgd,whichGO = go_results$`Term Id`[x])[[1]]),collapse = ", ")
      }))
      
      # removeModal()
      go_results
      
    })
    
    output$go_table <- renderReactable({
      
      req(go_results())
      data <- go_results()
      names(data) <- c("ID","Term","Reference","Gene_total","Enrichment","adj_pval","Genes")
      data_ <- data[,c("ID","Term","Gene_total","Enrichment","adj_pval","Genes")]
      reactable(
        data_,
        highlight = TRUE,
        bordered = TRUE,
        compact = TRUE,
        filterable = TRUE,
        showSortable = TRUE,
        defaultColDef = colDef(footer = function(values) {
          if (!is.numeric(values)) return()
          sparkline(values,type = "box",width = 100,height = 30)
        }),
        columns = list(
          Gene_total = colDef(
            name = "Total Genes"
            # ,
            # align = "left",cell = function(value) {
            # width <- paste0(value / max(data$Reference) * 100,"%")
            # bar_chart(value,width = width,fill = "#fc5185",background = "#e1e1e1")
            # }
          ),
          Genes = colDef(
            html = TRUE,
            cell =  function(value,index,name) {
              render.reactable.cell.with.tippy(text = value,tooltip = value)
          }),
          Enrichment = colDef(
            filterable = TRUE,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] >= filterValue
              })
            }")
          ),
          adj_pval = colDef(
            name = "Adj. p-value",
            filterable = TRUE,
            filterMethod = JS("function(rows, columnId, filterValue) {
              return rows.filter(function(row) {
                return row.values[columnId] >= filterValue
              })
            }")
          )
        )
      )
      
    })
    
    output$download_enrichment <- downloadHandler(
      
      filename = function() {
        'enrichment_results.txt'
      },
      content = function(file) {
        write.table(go_results(),file,sep = "\t",row.names = FALSE,col.names = TRUE)
      }
      
    )
    
}

# Run the application
shinyApp(ui = ui, server = server)
