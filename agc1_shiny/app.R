library(plotly)
library(DT)
library(tibble)
library(dplyr)
library(ggplot2)
library(shiny)

res <- readRDS("data_app/002_res.rds") |> 
  as.data.frame() |> 
  rownames_to_column("ensembl_id") |> 
  dplyr::rename(log2FC = log2FoldChange) 

ui <- fluidPage(
  
 fluidRow(
   column(2, 
          textInput("gene_vec", "Comma-separated gene symbols", placeholder = "Srebf1"),
          numericInput("logfc_th", "Log2Fold Change threshold", value = 1),
          numericInput("pth", "Padj threshold", value = 0.05, min = 0, max = 1)
          ),
   column(3, div(DTOutput("dge")), style = "font-size:80%"),
   column(6, plotlyOutput("vol"))
  )
)

server <- function(input, output, session) {
  
  # Create reactive expression that filters dataframe based on character vector input
  filtered_res <- reactive({
    if (input$gene_vec != "") {
      gene_vec_split <- strsplit(gsub(" ", "", input$gene_vec), ",")[[1]]
      res[res$symbol %in% gene_vec_split,]
    } else {
      return(res)
    }
  })
  output$dge <- renderDT(filtered_res() |> 
                           select(symbol, ensembl_id, log2FC, padj) |> 
                           relocate(symbol) |> 
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
  output$vol <- renderPlotly({
    
    lt <- input$logfc_th
    pt <- input$pth
    filtered_res() |> 
      filter(!is.na(padj)) |> 
      mutate(
        logp = -log10(padj),
        across(where(is.numeric), ~round(.x, 4)),
        color = case_when(
          abs(log2FC) > lt & logp > -log10(pt) ~ "p-value and log2FC",
          abs(log2FC) > lt & logp < -log10(pt) ~ "log2FC",
          abs(log2FC) < lt & logp > -log10(pt) ~ "value",
          abs(log2FC) < lt & logp < -log10(pt) ~ "NS"
        )) |> 
      ggplot(aes(log2FC, logp, color = color)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_vline(xintercept = c(-lt, lt), lty = 2, lwd = 0.5) +
      geom_hline(yintercept = pt, lty = 2, lwd = 0.5) +
      theme_light() +
      scale_y_continuous(limits = c(0, 50), oob = scales::squish) +
      # scale_x_continuous(limits = c(-5, 5), oob = scales::squish) +
      scale_color_manual(values = c("p-value and log2FC"= "red",
                                    "log2FC" = "#218c21",
                                    "value" = "#3f68e1",
                                    "NS" = "#4b4b4b"),
                         name = "",
                         na.translate = F)
    ggplotly() |> 
      layout(legend = list(orientation = "h", 
                           y = 1, x = 0.5, 
                           yanchor = "bottom", xanchor = "center"),
             font = list(family = "Arial"),
             xaxis = list(title = paste0("Log", "<sub>", "2", "</sub>", "FoldChange")),
             yaxis = list(title = paste0("-Log", "<sub>", "10", "</sub>", "<i>", "P", "</i>"))
             ) |> 
      toWebGL()
    
  })
}
shinyApp(ui, server)
