# Interactive version of the main results. Run from the repo root:
#   Rscript -e "shiny::runApp('shiny')"
# Needs shiny and ggplot2: install.packages(c("shiny", "ggplot2"))

library(shiny)
library(ggplot2)

perf <- read.csv("../results/reports/09_performance.csv")
prop <- read.csv("../results/reports/07_proportions.csv")
de_summary <- read.csv("../results/reports/08_de_summary.csv")

# shorter labels for the plots
perf$setting <- sub(" only", "", perf$setting)
prop$cell_type <- sub(", alpha-beta T cell", " T cell", prop$cell_type)
de_summary$cell_type <- sub(", alpha-beta T cell", " T cell", de_summary$cell_type)

ui <- fluidPage(
    titlePanel("Can a blood sample tell who has lupus?"),
    p("1.26M single cells from 162 lupus patients and 99 healthy donors.",
      "The pooled 0.99 AUROC is inflated: batch and ancestry alone already get 0.87.",
      "Batch 4 is the only balanced batch, and there the full model gets 0.92."),
    h4("Classifier AUROC, pooled and within batch"),
    plotOutput("auroc", height = "480px"),
    p("Note: the covariate model (batch + ancestry) has no bar for the single-batch",
      "settings. It cannot be trained when batch is constant."),
    fluidRow(
        column(6,
            h4("Cell-type composition"),
            p("Odds ratios for the SLE vs healthy mix of each cell type,",
              "95% CI included. Points right of 1 are expanded in lupus blood."),
            plotOutput("forest", height = "450px")
        ),
        column(6,
            h4("Hallmark pathways by cell type"),
            p("Interferon alpha and gamma response lead in every cell type (camera,",
              "pseudobulk limma-voom)."),
            imageOutput("pathways")
        )
    ),
    h4("Pseudobulk differential expression per cell type"),
    tableOutput("de_table")
)

server <- function(input, output) {
    output$auroc <- renderPlot({
        perf$setting <- factor(perf$setting)
        ggplot(perf, aes(x = setting, y = AUROC, fill = model)) +
            geom_col(position = position_dodge(0.8), width = 0.7) +
            geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                          position = position_dodge(0.8), width = 0.3, colour = "grey40") +
            scale_y_continuous(limits = c(0, 1.05)) +
            labs(x = NULL, y = "AUROC", fill = NULL) +
            theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  legend.position = "top")
    })

    output$forest <- renderPlot({
        ggplot(prop, aes(x = odds_ratio, y = reorder(cell_type, odds_ratio))) +
            geom_pointrange(aes(xmin = ci_low, xmax = ci_high), colour = "#4C72B0") +
            geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
            scale_x_log10() +
            labs(x = "odds ratio (SLE vs healthy), log scale", y = NULL) +
            theme_minimal(base_size = 12)
    })

    output$pathways <- renderImage({
        list(src = "../results/figures/08_pathways_heatmap.png",
             width = "100%", contentType = "image/png")
    }, deleteFile = FALSE)

    output$de_table <- renderTable({
        d <- de_summary[, c("cell_type", "de_up", "de_down", "deseq2_rho",
                            "pct_limma_de_also_deseq2", "top_up", "top_down")]
        names(d) <- c("cell type", "DE up", "DE down", "DESeq2 rho",
                      "limma also DESeq2", "top up genes", "top down genes")
        d$`DESeq2 rho` <- round(d$`DESeq2 rho`, 2)
        d$`limma also DESeq2` <- paste0(round(d$`limma also DESeq2`, 1), "%")
        d
    }, rownames = FALSE)
}

shinyApp(ui, server)
