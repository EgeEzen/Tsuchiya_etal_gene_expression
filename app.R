# renv::install("dplyr")
# renv::install("ggplot2")
# renv::install("ggthemr")
# renv::install("readxl")
# renv::install("shiny")
# renv::install("tidyr")
# renv::install("rsconnect")
# renv::install("ggtext")
# renv::restore()

#> shiny::runApp("~/Desktop/PhD/2024 Second Term/Shiny app")
#> rsconnect::deployApp("~/Desktop/PhD/2024 Second Term/Shiny app",appFileManifest = "manifest.txt")

library(shiny)
library(ggplot2)
library(ggthemr)
library(readxl)
library(tidyr)
library(dplyr)
library(ggtext)

#read the data for dot plot
results_list_oa <- readRDS("results_list_oa.rds")
results_list_ra <- readRDS("results_list_ra.rds")
gene_list <- readRDS("gene_list.rds")

#read the data for violin plot
violin_data_oa <- readRDS("violin_oa_data.rds")
violin_data_ra <- readRDS("violin_ra_data.rds")
oa_significance_lookup <- readRDS("oa_significance_lookup.rds")
ra_significance_lookup <- readRDS("ra_significance_lookup.rds")

# UI and server code follow below

# UI code
ui <- fluidPage(
  titlePanel("Gene Expression Tsuchiya"),
    tabsetPanel(
      tabPanel(
        title="Dotplot",
        sidebarLayout(
          sidebarPanel(
            selectizeInput("gene", "Select Gene Name:", choices = NULL, multiple = FALSE, 
                           options = list(placeholder = "Type to search gene...")),
            
            selectInput("condition_type", "Choose Condition:", choices = c("RA", "OA"), selected = "RA"),
            actionButton("update", "Draw Plot"),
            tags$hr(),
            numericInput("y_min", "Min Y-axis Limit", value = 0, step = 0.1),
            numericInput("y_max", "Max Y-axis Limit", value = 0, step = 0.1),
            tags$hr(),
            div(
              actionButton("toggleSettings", "Settings for Download", icon = icon("cogs")),
              conditionalPanel(
                condition = "input.toggleSettings % 2 == 1",  # Toggle visibility based on button click
                numericInput("plot_width", "Plot Width (inches)", value = 6.5, min = 3, max = 20, step = 0.5),
                numericInput("plot_height", "Plot Height (inches)", value = 4, min = 3, max = 20, step = 0.5)
              )
            ),
            downloadButton("downloadPlot", "Download Plot as PNG"),
            tags$hr(),
            downloadButton("downloadData", "Download Data as CSV") 
          ),
          mainPanel(
            plotOutput("dotplot")
          )
        )
      ),
      
      tabPanel(
        title = "Violin Plot",
        sidebarLayout(
          sidebarPanel(
            selectizeInput("violin_gene", "Select Gene Name:", choices = NULL, multiple = FALSE, 
                           options = list(placeholder = "Type to search gene...")),
            
            selectInput("violin_condition_type", "Choose Condition:", choices = c("RA", "OA"), selected = "RA"),
            actionButton("violin_update", "Draw Violin Plot"),
            tags$hr(),
            div(
              actionButton("toggleSettings", "Settings for Download", icon = icon("cogs")),
              conditionalPanel(
                condition = "input.toggleSettings % 2 == 1",  # Toggle visibility based on button click
                numericInput("violin_width", "Plot Width (inches)", value = 6.5, min = 3, max = 20, step = 0.5),
                numericInput("violin_height", "Plot Height (inches)", value = 4, min = 3, max = 20, step = 0.5)
              )
            ),
            downloadButton("downloadViolinPlot", "Download Violin Plot as PNG"),
            tags$hr(),
            downloadButton("downloadViolinData", "Download Data as CSV")
          ),
          mainPanel(
            plotOutput("violinPlot")
          )
        )
      )
    ),
  tags$hr(), 
  fluidRow(
    column(12,
           tags$footer(
             HTML("<p style='text-align: center; color: #666;'>
          This Shiny app visualizes gene expression data from
          <a href=\"https://pubmed.ncbi.nlm.nih.gov/33139312/\">Tsuchiya et. al.</a> 
          for RA and OA conditions across different stimulations ('all' stands for 8-mix). The publicly available 
          raw count data was analyzed with DESeq2. <br>
          Use the Dotplot tab to view expression changes. When the arguments for ylim and xlim are 0, the limits
          are determined automatically. <br>
          Use the Violin Plot tab to see the distribution of expression values. 
          When the change compared to non-stimulated (NS) is significant
          (padj < 0.05) the x-axis label for that stimulation will be bold. <br>
          For questions or feedback, you can contact me: ege.ezen@usz.ch. :)
        </p>
        <p style='text-align: center; font-size: 0.8em;'>
          2024 Ege Ezen / University of Zurich - Department of Rheumatology / Caroline Ospelt Lab <br>
          Last updated: December 2024.
        </p>"),
             style = "margin-top: 20px; padding: 10px; background-color: #f8f8f8;"
             
           )
    )
  )
)

# Server code
server <- function(input, output, session) {
  
  # --- DOTPLOT TAB ---
  
  # Initialize server-side selectize input
  updateSelectizeInput(session, "gene", choices = gene_list, server = TRUE)
  
  # Reactive expression to return the selected dataset (RA or OA)
  selected_results_list <- reactive({
    if (input$condition_type == "RA") {
      return(results_list_ra)
    } else {
      return(results_list_oa)
    }
  })

  observeEvent(input$update, {
    results_list <- selected_results_list()
    selected_gene <- input$gene
    
    # Initialize an empty dataframe for the dot plot
    dotplot_df <- data.frame(log2FoldChange = numeric(), padj = numeric(), condition = character())
    
    gather_gene_for_dotplot <- function(df, string) {
      append_to_dotplot_df <- df[df$genename %in% string, c("log2FoldChange", "padj")]
      return(append_to_dotplot_df)
    }
    
    for (name in names(results_list)) {
      df <- results_list[[name]]
      result <- gather_gene_for_dotplot(df, selected_gene)
      
      if (nrow(result) > 0) {
        result$condition <- name
        dotplot_df <- rbind(dotplot_df, result)
      }
    }
    
    if (nrow(dotplot_df) == 0) {
      output$dotplot <- renderPlot({
        ggplot() + 
          geom_text(aes(x = 1, y = 1, label = "Gene not found!"), size = 6, color = "red") +
          theme_void()
      })
    } else {
      dotplot_df$padj[is.na(dotplot_df$padj)] <- 0.999
      min_logFC <- min(dotplot_df$log2FoldChange, na.rm = TRUE)
      max_logFC <- max(dotplot_df$log2FoldChange, na.rm = TRUE)
      
      if (input$y_min == 0 && input$y_max == 0) {
        y_min <- (min_logFC) - 0.1
        y_max <- (max_logFC) + 0.1
      } else {
        y_min <- input$y_min
        y_max <- input$y_max
      }
      
      # Render the dot plot with ggthemr theme
      output$dotplot <- renderPlot({
        ggthemr('pale', layout= 'clear')  # Apply the ggthemr theme
        
        ggplot(dotplot_df, aes(x = condition, y = log2FoldChange, size = -log10(padj), 
                               color = ifelse(log2FoldChange > 0 & padj < 0.05, "logFC > 0", 
                                              ifelse(log2FoldChange < 0 & padj < 0.05, "logFC < 0", "Not Significant")))) +
          geom_point(alpha = 1) +
          scale_size_continuous(range = c(4, 7), breaks = c(0.3, 0.9, 3)) +
          labs(x = "Condition", y = "logFC", size = "-log10(padj)", title = paste("Gene Expression for", selected_gene, "in", input$condition_type)) +
          ylim(c(y_min, y_max)) +  # Set the y-axis limits dynamically or manually
          geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
          theme(
            legend.box.background = element_rect(color = "darkseagreen"),  # Add box around legend
            legend.box.margin = margin(5, 5, 5, 5), ## Set margin around legend box
            axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
          scale_color_manual(values = c("logFC > 0" = "#E97777", "logFC < 0" = "#4D869C", "Not Significant" = "grey")) +
          labs(color = "LogFC Colors", size = "-log10(padj)") +
          guides(size = guide_legend(override.aes = list(color = "grey")))  # Change the size and color of points in the size legend
      })
    }
  
    # Create a download handler for the data
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(selected_gene, "_dotplot_data.csv", sep = "")
      },
      content = function(file) {
        write.csv(dotplot_df, file, row.names = FALSE)
      }
    )
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(input$gene,"_" ,input$condition_type, "_dotplot.png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, bg = "transparent", 
             width = input$plot_width, height = input$plot_height)
    }
  )
  
  
  # --- VIOLIN PLOT TAB ---
  
  # Update selectize input dynamically with gene list for Violin Plot
  updateSelectizeInput(session, "violin_gene", choices = gene_list, server = TRUE)
  
  # Expression for Violin Plot
  observeEvent(input$violin_update, {
    gene_of_interest <- input$violin_gene
    condition_type <- input$violin_condition_type
    
    # Filter the data for the selected gene and condition type
    selected_violin_data <- if (condition_type == "OA") {
      violin_data_oa
    } else {
      violin_data_ra
    }
    
    selected_violin_lookup <- if (condition_type == "OA") {
      oa_significance_lookup
    } else {
      ra_significance_lookup
    }
    
    condition_vector <- c("NS", "all", "IFNa", "IFNg", "IL17", "IL18", "IL1b", "IL6", "TGFb", "TNFa")
    
    significant_data <- selected_violin_lookup %>%
      filter(genename == gene_of_interest) %>%
      filter(significance == "S") %>%
      select(gene = genename, stimulation, significance)
    
    significant_stimulation <- significant_data$stimulation
    
    reshaped_filtered <- selected_violin_data %>%
      filter(gene == gene_of_interest) %>%
      mutate(
        stimulation = factor(stimulation, levels = condition_vector),
        StyledStimulation = ifelse(stimulation %in% significant_stimulation, 
                                   paste0("<b>", stimulation, "</b>"), 
                                   as.character(stimulation))
      ) %>%
      arrange(stimulation) %>% #so that it doesn't lose the factor
      mutate(StyledStimulation = factor(StyledStimulation, levels = unique(StyledStimulation), ordered = TRUE))
    
    # Define condition colors
    condition_colors <- c(
      "NS" = "#FAD02E",
      "All" = "#F28D35",
      "IFNa" = "#F77F00",
      "IFNg" = "#D4A5A5",
      "IL17" = "#96B1E6",
      "IL18" = "#F4C2C2",
      "IL1b" = "#D1F0C9",
      "IL6" = "#C5B9D9",
      "TGFb" = "#9ECAE1",
      "TNFa" = "#FFD1DC"
    )
    
    # Generate the violin plot
    if (nrow(reshaped_filtered)== 0){
      showNotification("Not enough counts to draw a violin plot.", type = "error")
      output$violinPlot <- renderPlot({
        ggplot() + 
          geom_text(aes(x = 1, y = 1, label = "No Data to Plot"), size = 6, color = "red") +
          theme_void()  # Creates an empty plot with a message
      })
    }
    else {
      output$violinPlot <- renderPlot({
      ggthemr('pale', layout = "clean")  # Apply theme
      
      ggplot(reshaped_filtered, aes(x = StyledStimulation, y = counts, fill = stimulation)) +
        geom_violin(trim = TRUE) +  # Violin plot
        geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color = "darkgrey") +  # Jittered dots
        labs(title = paste("Gene Expression for", gene_of_interest, "in", condition_type), 
             x = "Stimulation", y = "Normalized Counts") +
        scale_fill_manual(values = condition_colors) +  # Apply custom colors
        theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1), 
              legend.position = "None")  # Rotate x-axis labels and remove legend
    })
    }
    # Create a download handler for the violin plot data
    output$downloadViolinData <- downloadHandler(
      filename = function() {
        paste(gene_of_interest, "_violin_data.csv", sep = "")
      },
      content = function(file) {
        write.csv(reshaped_filtered, file, row.names = FALSE)
      }
    )
    
    # Create a download handler for the violin plot as PNG
    output$downloadViolinPlot <- downloadHandler(
      filename = function() {
        paste(input$violin_gene, "_violin_plot.png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png", dpi = 300, bg = "transparent", 
               width = input$violin_width, height = input$violin_height)
      }
    )
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)

