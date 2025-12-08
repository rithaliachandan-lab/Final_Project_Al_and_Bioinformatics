# ================================
# TCGA LUAD SHINY DASHBOARD
# Gene Expression + Correlation + Heatmap + PCA
# ================================

library(shiny)

# -------------------------------
# Load processed TCGA data
# -------------------------------
ss_rna <- readRDS("ss_rna_for_shiny.rds")

# -------------------------------
# Ensure log2 columns exist
# -------------------------------
if (!"log2_TPM_NKX2_1" %in% colnames(ss_rna)) {
  ss_rna$log2_TPM_NKX2_1 <- log2(ss_rna$TPM_NKX2_1 + 1)
}
if ("TPM_EGFR" %in% colnames(ss_rna) && !"log2_TPM_EGFR" %in% colnames(ss_rna)) {
  ss_rna$log2_TPM_EGFR <- log2(ss_rna$TPM_EGFR + 1)
}
if ("TPM_TP53" %in% colnames(ss_rna) && !"log2_TPM_TP53" %in% colnames(ss_rna)) {
  ss_rna$log2_TPM_TP53 <- log2(ss_rna$TPM_TP53 + 1)
}

# -------------------------------
# Gene → column mapping
# -------------------------------
gene_map <- c(
  "NKX2-1" = "log2_TPM_NKX2_1",
  "EGFR"   = "log2_TPM_EGFR",
  "TP53"   = "log2_TPM_TP53"
)

gene_choices <- names(gene_map)

# ===============================
# UI
# ===============================
ui <- fluidPage(
  titlePanel("TCGA LUAD Bioinformatics Dashboard"),
  
  tabsetPanel(
    
    # -------- TAB 1: Expression --------
    tabPanel("Expression",
             sidebarLayout(
               sidebarPanel(
                 selectInput("expr_gene", "Select gene:",
                             choices = gene_choices,
                             selected = "NKX2-1"),
                 
                 checkboxGroupInput("expr_types",
                                    "Sample types:",
                                    choices = c("Primary Tumor", "Solid Tissue Normal"),
                                    selected = c("Primary Tumor", "Solid Tissue Normal"))
               ),
               mainPanel(
                 plotOutput("expr_boxplot")
               )
             )
    ),
    
    # -------- TAB 2: Correlation --------
    tabPanel("Correlation (Tumor only)",
             sidebarLayout(
               sidebarPanel(
                 selectInput("cor_x_gene", "X-axis gene:",
                             choices = gene_choices,
                             selected = "NKX2-1"),
                 selectInput("cor_y_gene", "Y-axis gene:",
                             choices = gene_choices,
                             selected = "EGFR")
               ),
               mainPanel(
                 plotOutput("cor_scatter"),
                 verbatimTextOutput("cor_stats")
               )
             )
    ),
    
    # -------- TAB 3: Heatmap --------
    tabPanel("Heatmap (Tumor only)",
             mainPanel(
               plotOutput("heatmap_plot", height = "500px")
             )
    ),
    
    # -------- TAB 4: PCA --------
    tabPanel("PCA (Tumor only)",
             sidebarLayout(
               sidebarPanel(
                 checkboxGroupInput("pca_genes",
                                    "Genes to include:",
                                    choices = gene_choices,
                                    selected = gene_choices)
               ),
               mainPanel(
                 plotOutput("pca_plot")
               )
             )
    )
  )
)

# ===============================
# SERVER
# ===============================
server <- function(input, output, session) {
  
  # Tumor-only data
  tumor_data <- reactive({
    subset(ss_rna, Sample.Type == "Primary Tumor")
  })
  
  # -------- EXPRESSION TAB --------
  output$expr_boxplot <- renderPlot({
    d <- ss_rna[ss_rna$Sample.Type %in% input$expr_types, ]
    colname <- gene_map[[input$expr_gene]]
    
    if (!colname %in% colnames(d) || all(is.na(d[[colname]]))) {
      plot.new()
      title("No data available for this gene")
      return()
    }
    
    boxplot(
      d[[colname]] ~ d$Sample.Type,
      xlab  = "Sample Type",
      ylab  = "log2(TPM + 1)",
      main  = paste(input$expr_gene, "Expression in LUAD"),
      col   = c("tomato", "skyblue")[seq_along(unique(d$Sample.Type))]
    )
  })
  
  # -------- CORRELATION TAB --------
  output$cor_scatter <- renderPlot({
    td <- tumor_data()
    x_col <- gene_map[[input$cor_x_gene]]
    y_col <- gene_map[[input$cor_y_gene]]
    
    if (x_col == y_col) {
      plot.new()
      title("Select two different genes")
      return()
    }
    
    td <- td[!is.na(td[[x_col]]) & !is.na(td[[y_col]]), ]
    if (nrow(td) < 3) {
      plot.new()
      title("Not enough samples")
      return()
    }
    
    plot(td[[x_col]],
         td[[y_col]],
         xlab = paste(input$cor_x_gene, "(log2 TPM+1)"),
         ylab = paste(input$cor_y_gene, "(log2 TPM+1)"),
         main = paste(input$cor_x_gene, "vs", input$cor_y_gene),
         pch  = 19)
  })
  
  output$cor_stats <- renderPrint({
    td <- tumor_data()
    x_col <- gene_map[[input$cor_x_gene]]
    y_col <- gene_map[[input$cor_y_gene]]
    
    td <- td[!is.na(td[[x_col]]) & !is.na(td[[y_col]]), ]
    if (nrow(td) < 3) {
      cat("Not enough tumor samples for correlation.")
      return()
    }
    
    print(cor.test(td[[x_col]], td[[y_col]], method = "pearson"))
  })
  
  # -------- HEATMAP TAB --------
  output$heatmap_plot <- renderPlot({
    td <- tumor_data()
    
    cols_needed <- c("log2_TPM_NKX2_1", "log2_TPM_EGFR", "log2_TPM_TP53")
    td <- td[complete.cases(td[, cols_needed]), ]
    
    if (nrow(td) < 2) {
      plot.new()
      title("Not enough tumor samples")
      return()
    }
    
    mat <- as.matrix(td[, cols_needed])
    colnames(mat) <- c("NKX2-1", "EGFR", "TP53")
    rownames(mat) <- NULL
    
    my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
    
    heatmap(mat,
            Rowv = NA,
            Colv = NA,
            scale = "column",
            col = my_colors,
            labRow = FALSE,
            main = "Heatmap of LUAD Biomarkers",
            xlab = "Genes",
            ylab = "Tumor Samples")
  })
  
  # -------- PCA TAB --------
  output$pca_plot <- renderPlot({
    td <- tumor_data()
    chosen_genes <- input$pca_genes
    
    if (length(chosen_genes) < 2) {
      plot.new()
      title("Select at least 2 genes")
      return()
    }
    
    cols <- gene_map[chosen_genes]
    td <- td[complete.cases(td[, cols]), ]
    
    if (nrow(td) < 3) {
      plot.new()
      title("Not enough samples for PCA")
      return()
    }
    
    pca_mat <- as.matrix(td[, cols])
    pca_res <- prcomp(pca_mat, scale. = TRUE)
    
    pc1_var <- round(summary(pca_res)$importance[2,1] * 100, 2)
    pc2_var <- round(summary(pca_res)$importance[2,2] * 100, 2)
    
    plot(pca_res$x[,1],
         pca_res$x[,2],
         pch = 19,
         col = "blue",
         xlab = paste0("PC1 (", pc1_var, "%)"),
         ylab = paste0("PC2 (", pc2_var, "%)"),
         main = "PCA of LUAD Tumor Samples")
  })
}

# ===============================
# RUN APP
# ===============================
shinyApp(ui, server)
