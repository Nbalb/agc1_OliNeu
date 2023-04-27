library(ggpointdensity)
library(tidyr)
library(dqshiny)
library(plotly)
library(DT)
library(dplyr)
library(ggplot2)
library(shiny)

# Load data
res <- readRDS("data/002_res.rds") |> 
  select(symbol, ensembl_id, log2FC, padj) |> 
  relocate(symbol)
msdb <- readRDS("data/msigdbr_reduced_db.rds") |> 
  group_by(gs_name) |> 
  mutate(n = n()) |> 
  ungroup() |> 
  filter(n > 20)
unique_path <- unique(msdb$gs_name)

tpms <- readRDS("data/002_tpms_tbl.rds") |> 
  mutate(genotype = gsub(".$", "", sample),
         across(where(is.numeric), ~round(.x, 3))) |> 
  rename(ensembl_id = gene_id) |>
  relocate(symbol, ensembl_id, log2FC, padj)

atac_res <- read.csv("data/ATAC_de.csv") |> 
  as_tibble() |> 
  select(gene_symbol, log2FoldChange, padj, geneChr, ensembl_gene) |> 
  rename(log2FC_peaks = log2FoldChange) |> 
  relocate(gene_symbol, ensembl_gene)
  
# Define UI
ui <- fluidPage(
  titlePanel(title = div(img(src = "Picture5.png", height = "100%", width = "100%"), 
                         style = "margin-left:-20px;margin-right:-20px;margin-top:-20px;"),
             windowTitle = "Agc1 DGE Explorer"),
 fluidRow(
   sidebarLayout(
     sidebarPanel(
       width = 2,
       textInput("gene_vec", "Comma-separated gene symbols", placeholder = "Srebf1, Srebf2, Slc25a12"),
       autocomplete_input("pathway", "Search MSigDB pathways", placeholder = "GOBP, KEGG, REACTOME", 
                          options = unique_path, contains = T),
       numericInput("logfc_th", "Log2Fold Change threshold", value = 1),
       numericInput("pth", "Padj threshold", value = 0.05, min = 0, max = 1)
       ),
   mainPanel(
     width = 10,
     fluidRow(
       tabsetPanel(
         type = "tabs",
         tabPanel("Volcano", 
                  column(4, div(DTOutput("dge"), 
                                # style = "font-size:90%"
                                )),
                  column(7, offset = 1, div(plotlyOutput("vol"), 
                                            style = "margin-right:30px"))
                  ),
         tabPanel("Boxplot", 
                  column(5, div(DTOutput("table_tpms"), 
                                # style = "font-size:90%"
                                )),
                  column(6, offset = 1, div(plotlyOutput("boxplot"),
                                            style = "margin-right:30px"
                                            ))
                  ),
         tabPanel("ATAC vs RNA", 
                  column(5, div(DTOutput("atac_table"), 
                                # style = "font-size:90%"
                  )),
                  column(6, offset = 1, div(plotlyOutput("atac_point"),
                                            style = "margin-right:30px"
                  ))
         )
         )
       )
     )
   )
   )
 )

server <- function(input, output, session) {
  pathway_ids <- reactive({
    msdb[msdb$gs_name %in% input$pathway,]$ensembl_gene
  })
  gene_vec_split <- reactive({
    strsplit(gsub(" ", "", input$gene_vec), ",")[[1]]
    })
  
    # Create reactive expression that filters dataframe based on character vector input
    filtered_res <- reactive({
      
      # Select genes based on the pathways they belong to
      if (input$pathway != "") {
        temp_res <- res[res$ensembl_id %in% pathway_ids(),]
      } else {
        temp_res <- res
      }
      
      # Select genes from text input
      if (input$gene_vec != "") {
        gene_vec_split <- temp_res[temp_res$symbol %in% gene_vec_split(),]
      } else {
        return(temp_res)
      }
    })
    
    # DataTable ----
    output$dge <- renderDT(filtered_res() |> 
                             mutate(across(where(is.numeric), ~round(.x, 3))), 
                           rownames = FALSE,
                           escape = FALSE,
                           options = list(
                             autoWidth = T,
                             order = list(list(3, "asc"), list(2, "desc")),
                             columnDefs = list(list(
                               targets = 1,
                               render = JS(
                                 "function(data, type, row, meta) {",
                                 "return type === 'display' && data.length > 6 ?",
                                 "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
                                 "}")
                             ))), 
                           callback = JS('table.page(3).draw(false);')
    )
    
    # Select genes from table click selection
    plot_res <- reactive({
      selected_rows <- input$dge_rows_selected
      if (!is.null(selected_rows)) {
        filtered_res()[selected_rows, , drop = FALSE]
      } else {
        return(filtered_res())
      }  
    })
    
    # Volcano Plot ----
    output$vol <- renderPlotly({
      
      lt <- input$logfc_th
      pt <- input$pth
      plot_res() |> 
        mutate(
          logp = -log10(padj),
          across(where(is.numeric), ~round(.x, 4)),
          color = case_when(
            abs(log2FC) > lt & logp > -log10(pt) ~ "p-value and log2FC",
            abs(log2FC) > lt & logp < -log10(pt) ~ "log2FC",
            abs(log2FC) < lt & logp > -log10(pt) ~ "p-value",
            (abs(log2FC) < lt & logp < -log10(pt)) | is.na(logp) ~ "NS"
          )) |> 
        ggplot(aes(log2FC, logp, color = color, 
                   text = paste0(symbol, "<br>",
                                "Log<sub>2</sub>FC: ", log2FC, "<br>",
                                "<i>P</i> adjusted: ", padj, "<br>"))) +
        geom_point(alpha = 0.5, size = 1.5) +
        geom_vline(xintercept = c(-lt, lt), lty = 2, lwd = 0.5) +
        geom_hline(yintercept = pt, lty = 2, lwd = 0.5) +
        theme_light() +
        scale_y_continuous(limits = c(0, 50), oob = scales::squish) +
        # scale_x_continuous(limits = c(-5, 5), oob = scales::squish) +
        scale_color_manual(values = c("p-value and log2FC"= "red",
                                      "log2FC" = "#218c21",
                                      "p-value" = "#3f68e1",
                                      "NS" = "#4b4b4b"),
                           name = "",
                           na.translate = F)
      ggplotly(tooltip = c("text"),
               height = 500) |>  
        config(toImageButtonOptions = list(format = "png", 
                                           scale = 6),
               displaylogo = F) |> 
        layout(legend = list(orientation = "h", 
                             y = 1, x = 0.5, 
                             yanchor = "bottom", xanchor = "center"),
               font = list(family = "Arial"),
               xaxis = list(title = paste0("Log", "<sub>", "2", "</sub>", "FoldChange")),
               yaxis = list(title = paste0("-Log", "<sub>", "10", "</sub>", "<i>", "P", "</i>"))
               ) |> 
        toWebGL()
    
  })
    
    # DataTable TPMS
    # Select genes based on the pathways they belong to
    filtered_tpms <- reactive({
      if (input$pathway != "") {
        temp_tpms <- tpms[tpms$ensembl_id %in% pathway_ids(),]
      } else {
        temp_tpms <- tpms
      }
      
      # Select genes from text input
      if (input$gene_vec != "") {
        temp_tpms <- temp_tpms[temp_tpms$symbol %in% gene_vec_split(),]
      }
      
      temp_tpms <- temp_tpms |> 
        select(-genotype) |> 
        pivot_wider(id_cols = everything(), 
                    names_from = sample,
                    values_from = tpms) 
      
      return(temp_tpms)
      
    })
   
    output$table_tpms <- renderDT(filtered_tpms(),
                                  rownames = FALSE,
                                  selection = "multiple", 
                                  escape = FALSE,
                                  options = list(
                                    autoWidth = T,
                                    order = list(list(3, "asc"), list(2, "desc")),
                                    columnDefs = list(list(
                                      targets = 1,
                                      render = JS(
                                        "function(data, type, row, meta) {",
                                        "return type === 'display' && data.length > 6 ?",
                                        "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
                                        "}")
                                    ))), 
                                  callback = JS('table.page(3).draw(false);')) 
    # Boxplot ----
    selection_id <- reactive({
      filtered_tpms()[input$table_tpms_rows_selected, ] |> 
        pull(ensembl_id) |> 
        as.character()
      })
    
    filt_long_tpms <- reactive({
      if(!is.null(input$table_tpms_rows_selected)){
        temp_long_tpms <- tpms[tpms$ensembl_id %in% selection_id(), ]
      } else if(input$pathway != "" | input$gene_vec != ""){
        temp_long_tpms <- tpms[tpms$ensembl_id %in% filtered_tpms()$ensembl_id, ]
        if(length(filtered_tpms()$ensembl_id) > 16){
          top16 <- filtered_tpms() |> 
            slice_max(abs(log2FC), n = 16) |> 
            pull(ensembl_id)
          temp_long_tpms <- temp_long_tpms[temp_long_tpms$ensembl_id %in% top16, ]
        }
      } else {
        temp_long_tpms <- tpms[tpms$symbol == "Slc25a12", ]
      }
      
      
      
      temp_long_tpms <- temp_long_tpms |> 
        arrange(desc(log2FC))
      temp_long_tpms |> 
        mutate(genotype = factor(genotype, 
                                 levels = c("wt", "kd"),
                                 labels = c("WT", "KD")),
               symbol = factor(symbol, levels = unique(temp_long_tpms$symbol)))
      })
    
    output$boxplot <- renderPlotly({

      filt_long_tpms() |>
        ggplot(aes(genotype, tpms, fill = genotype)) +
        geom_boxplot(alpha = 0.5, color = "black") +
        facet_wrap(~symbol, ncol = 4, scales = "free") +
        labs(x = "Genotype",
             y = "TPMS") +
        scale_fill_discrete(name = "Genotype") + 
        theme_bw() +
        theme(panel.spacing.x = unit(1, "lines"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"), 
              panel.border = element_rect(color = "black", fill = NA),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) 
      ggplotly(height = 500) |> 
        config(toImageButtonOptions = list(format = "png", 
                                           scale = 6),
               displaylogo = F) |> 
        layout(font = list(family = "Arial"),
               legend = list(x = 100,
                             y = 0.5))
      
    })

    # DT ATAC data
    filtered_atac <- reactive({
      
      # Select genes based on the pathways they belong to
      if (input$pathway != "") {
        temp_atac_res <- atac_res[atac_res$ensembl_id %in% pathway_ids(),]
      } else {
        temp_atac_res <- atac_res
      }
      
      # Select genes from text input
      if (input$gene_vec != "") {
        gene_vec_split <- temp_atac_res[temp_atac_res$symbol %in% gene_vec_split(),]
      } else {
        return(temp_atac_res)
      }
    })
    
    output$atac_table <- renderDT(filtered_atac() |> 
                             mutate(across(where(is.numeric), ~round(.x, 3))), 
                           rownames = FALSE,
                           escape = FALSE,
                           options = list(
                             autoWidth = T,
                             order = list(list(3, "asc"), list(2, "desc")),
                             columnDefs = list(list(
                               targets = 1,
                               render = JS(
                                 "function(data, type, row, meta) {",
                                 "return type === 'display' && data.length > 6 ?",
                                 "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
                                 "}")
                             ))), 
                           callback = JS('table.page(3).draw(false);')
    )
    
    # Select genes from table click selection
    plot_atac <- reactive({
      selected_rows <- input$atac_table_rows_selected
      if (!is.null(selected_rows)) {
        filtered_atac()[selected_rows, , drop = FALSE]
      } else {
        return(filtered_atac())
      }  
    })
    
    # ATAC plots ----
    output$atac_point <- renderPlotly({
      
      plot_atac() |> 
        inner_join(res |> filter(padj < 0.05),
                   by = c("gene_symbol" = "symbol"), suffix = c("_ATAC", "_RNA")) |> 
        mutate(across(where(is.numeric), ~round(.x, 3))) |> 
        rename(log2FC_RNA = log2FC, log2FC_ATAC = log2FC_peaks, symbol = gene_symbol) |> 
        ggplot(aes(log2FC_RNA, log2FC_ATAC, 
                   text = paste0(symbol, "<br>",
                                 "Log<sub>2</sub>FC RNA: ", log2FC_RNA, "<br>",
                                 "<i>P</i> adjusted RNA: ", padj_RNA, "<br>",
                                 "Log<sub>2</sub>FC ATAC: ", log2FC_ATAC, "<br>",
                                 "<i>P</i> adjusted ATAC: ", padj_ATAC, "<br>"))) +
        geom_point(alpha = .7, color = "black") + 
        labs(title = "ATAC differential accessibility vs RNA differential gene expression\n
             Significant results only",
             x = "RNA-Seq log2FoldChange",
             y = "ATAC-Seq log2FoldChange") + 
        theme_bw()  
      ggplotly(tooltip = "text", 
               height = 500) |> 
        config(toImageButtonOptions = list(format = "png", 
                                           scale = 6),
               displaylogo = F) |> 
        layout(font = list(family = "Arial"),
               legend = list(x = 100,
                             y = 0.5),
               title = list(text = paste0('ATAC differential accessibility vs RNA differential gene expression',
                                    '<br>',
                                    '<sup>',
                                    'Significant results only',
                                    '</sup>'))) |> 
        toWebGL()
    })
    
}
shinyApp(ui, server)
