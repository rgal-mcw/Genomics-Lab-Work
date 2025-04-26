#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(tidyverse)
library(shinyjs)

# FILE INFO - MANUAL INPUT AS OF RN
file = "~/Desktop/Project/maple/svi-0101_UDD/mcw_svi_0101_UDD_vcf_calls.tsv"
hets = read.table("~/Desktop/Project/maple/svi-0101_UDD/compound_hets.tsv", sep="\t", header=T, fill=T, comment.char="", quote="", stringsAsFactors = F)
sample_name = paste((strsplit(basename(file), "_")[[1]])[1:3], collapse = "_")
vcf = read.table(file, sep="\t", header=TRUE, fill=TRUE, comment.char="", quote="", stringsAsFactors = F)

ui = page_fluid(
  useShinyjs(), #Init shinyjs
  navset_tab(
    
    # Header content to include custom CSS and the centered title
    header = tagList(
      tags$head(
        tags$style(HTML("
          .navbar-custom-title {
            position: absolute;
            left: 50%;
            transform: translateX(-50%);
            font-weight: bold;
            font-size: 30px;
            top: 0px; /* Adjust this value to control vertical alignment */
            z-index: 1000;
          }
          .top-right-image {
            position: absolute;
            top: 10px;
            right: 10px;
            z-index: 1001;
            width: 100px; /* Adjust the size as needed */
          }
          .dataTable td {
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 300px;
            border-right: 1px solid #ddd;  /* Vertical borders between columns */
          }
          .dataTable th {
            border-right: 1px solid #ddd;  /* Vertical borders between headers */
          }
          .dataTable td:hover {
            cursor: pointer;
          }
          /* Global font styling */
          body {
            font-family: 'Franklin Gothic Book', sans-serif;
            font-size: 16px;
          }
          h1, h2, h3, h4, h5, h6 {
            font-family: 'Franklin Gothic', sans-serif;
          }
        "))
      ),
      
      # Centered title
      tags$div(
        paste(sample_name, "Interactive VCF"),
        class = "navbar-custom-title"
      ),
      # Image at the top right corner
      tags$img(
        src = "16bit_toucan_art.png",
        class = "top-right-image"
      )
    ),
      
    # Page 1
    nav_panel("Table", fluidPage(
        
              # Center the Show/Hide Options button
              div(
                style = "text-align: left; margin-top: 10px; margin-bottom: 10px",  # Center the button
                actionButton("toggle_btn", "Show/Hide Table Filters")
              ),
              
              # The input panel (this will be toggled)
              hidden(div(
                id = "input_panel",
                fluidRow(
                  column(4, 
                         selectInput(
                           inputId = "selected_vars", 
                           label = "Select variables to display:", 
                           # Dynamically reorder the choices so the selected columns come first (must do it both lists in choices)
                           choices = c(intersect(c("CHROM", "POS", "ID", "GT", "GENEINFO", "CLNSIG", "MC", "CLNDN", 
                                                   "CLNSIGCONF", "MC_details", "IMPACT_details", "GENE_details", "OMIM_INFO", "AF_EXAC", "AF_TGP", "AF_ESP"), 
                                                 colnames(vcf)), 
                                       setdiff(colnames(vcf), 
                                               c("CHROM", "POS", "ID", "GT", "GENEINFO", "CLNSIG", "MC", "CLNDN", 
                                                 "CLNSIGCONF", "MC_details", "IMPACT_details", "GENE_details", "OMIM_INFO","AF_EXAC", "AF_TGP", "AF_ESP"))
                           ),
                           
                           # Keep the default selected columns
                           selected = c("CHROM", "POS", "ID", "GT", "GENEINFO", "CLNSIG", "MC", "CLNDN", 
                                        "CLNSIGCONF", "MC_details", "IMPACT_details", "GENE_details", "OMIM_INFO", "AF_EXAC", "AF_TGP", "AF_ESP"),                            multiple = TRUE
                         )
                  ),
                  column(4, 
                         selectizeInput(
                           inputId = "filter_columns", 
                           label = "Select columns to filter non-empty rows:", 
                           choices = colnames(vcf),  
                           selected = c("CHROM", "POS"),  
                           multiple = TRUE  
                         ),
                         
                         # Checkbox for case-sensitive filtering (stacked under selectizeInput)
                         checkboxInput(
                           inputId = "case_sensitive", 
                           label = "Case-sensitive filtering", 
                           value = TRUE  
                         )
                  )
                )
              )),
              
              # Display the table
              DTOutput("data_table"),
        
              # Center the buttons and the generated links
              div(
                style = "display: flex; justify-content: center; flex-direction: column; align-items: center; margin-top: 20px;",  # Flexbox for vertical alignment
                
                # Buttons for Prepare NCBI and Prepare OMIM
                div(
                  style = "display: flex; gap: 10px;",  # Horizontal alignment of buttons with spacing
                  actionButton("go_to_url", "Prepare NCBI SNP Page"),
                  actionButton("go_to_omim", "Prepare OMIM Search Page")
                ),
                
                # Generated links for NCBI and OMIM (directly beneath the buttons)
                div(
                  style = "margin-top: 10px;",  # Add space between the buttons and links
                  uiOutput("snp_link"),
                  uiOutput("omim_link")
                )
              ),
              
              # Show the selected cell value
              div(
                style = "white-space: normal; word-wrap: break-word;",  # CSS to enable text wrapping
                textOutput("selected_cell")
              )
          )),
  
    # Page 2 - Compound Hets
    # Page 2 - Compound Hets
    nav_panel("Compound Hets", fluidPage(
      DTOutput("hets_table")
    )),
    
    # Page 3
    nav_panel("C", "Page C content"), 
    nav_menu( 
      "Other links", 
      nav_panel("D", "Panel D content"), 
      "----", 
      "Description:", 
      nav_item( 
        a("Shiny", href = "https://shiny.posit.co", target = "_blank") 
      ), 
    ), 
  ), 
  id = "tab" 
)

server <- function(input, output) {
  
        # Toggle the visibility of the input panel when the button is clicked
        observeEvent(input$toggle_btn, {
          toggle("input_panel")
        })
  
        # Create a reactive expression to filter out rows with empty cells in the selected columns
        filtered_data <- reactive({
          req(input$selected_vars, input$filter_columns)  # Ensure inputs are available
          
          # Get the subset of the dataframe with the selected columns
          filtered <- vcf[, input$selected_vars, drop = FALSE]
          
          # Filter out rows where any of the selected columns have empty or NA values
          for (column in input$filter_columns) {
            column_to_filter <- filtered[[column]]
            filtered <- filtered[!is.na(column_to_filter) & column_to_filter != "", ]
          }
          
          filtered  # Return the filtered dataframe
        })
        
        # Render the datatable with single cell selection and column-specific filtering
        output$data_table <- renderDT({
          datatable(
            filtered_data(),
            rownames = FALSE,          # Remove the index column
            extensions = 'ColReorder',
            filter = "top",            # Enable column-specific filtering
            selection = list(mode = "single", target = "cell"),  # Single-cell selection
            options = list(
              pageLength = 10,
              scrollX = TRUE,           # Enable horizontal scrolling
              colReorder = TRUE,        # Enable column reordering
              autoWidth = TRUE,         # Automatically adjust column widths
              #dom = 'tipfl', # Move "Show Entries" & Search to bottom
              search = list(
                caseInsensitive = !input$case_sensitive  # Toggle case-sensitive/insensitive filtering
              )
            )
          )
        })
        
        # Display the selected cell value
        output$selected_cell <- renderText({
          selected <- input$data_table_cells_selected
          if (length(selected)) {
            row <- selected[1]
            col <- selected[2] + 1  # Adjust for 1-based indexing
            cell_value <- filtered_data()[row, col]
            paste(cell_value)
          } else {
            "No cell selected"
          }
        })
        
        output$hets_table <- renderDT({
          datatable(
            hets,
            rownames = FALSE,
            extensions = 'ColReorder',
            filter = "top",
            selection = list(mode = "single", target = "cell"),
            options = list(
              pageLength = 10,
              scrollX = TRUE,
              colReorder = TRUE,
              autoWidth = TRUE
            )
          )
        })
        
        # Observe button click to update the link for NCBI SNP page
        observeEvent(input$go_to_url, {
          selected <- input$data_table_cells_selected
          if (length(selected)) {
            row <- selected[1]  # Get selected row
            col <- selected[2]  # Get selected column
            
            # Check if the selected cell is in the ID column (assuming rsID is in column "ID")
            if (colnames(filtered_data())[col + 1] == "ID") {
              rsID <- filtered_data()[row, col + 1]  # Adjust for 1-based indexing
              url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", rsID)  # Construct the full URL
              
              # Display the link for the constructed URL
              output$snp_link <- renderUI({
                tags$a(href = url, target = "_blank", "Click here to open NCBI SNP Page", style = "font-weight:bold;")
              })
            } else {
              showNotification("Please select an rsID cell (from the ID column).", type = "warning")
              output$snp_link <- renderUI(NULL)  # Hide the link
            }
          } else {
            showNotification("No cell selected!", type = "error")
            output$snp_link <- renderUI(NULL)  # Hide the link
          }
        })
        
        # Observe button click to update the link for OMIM gene page
        observeEvent(input$go_to_omim, {
          selected <- input$data_table_cells_selected
          if (length(selected)) {
            row <- selected[1]  # Get selected row
            col <- selected[2]  # Get selected column
            
            # Check if the selected cell is in the GENEINFO column (extract gene name)
            if (colnames(filtered_data())[col + 1] == "GENEINFO") {
              gene_info <- filtered_data()[row, col + 1]  # Get the GENEINFO value
              gene_name <- strsplit(gene_info, ":")[[1]][1]  # Extract text before the first ':'
              url <- paste0("https://www.omim.org/search?index=entry&start=1&limit=10&sort=score+desc%2C+prefix_sort+desc&search=", gene_name)  # Construct the OMIM URL
              
              # Display the link for the constructed OMIM URL
              output$omim_link <- renderUI({
                tags$a(href = url, target = "_blank", paste("Click here to open OMIM Search for", gene_name), style = "font-weight:bold;")
              })
            } else {
              showNotification("Please select a cell from the GENEINFO column.", type = "warning")
              output$omim_link <- renderUI(NULL)  # Hide the link
            }
          } else {
            showNotification("No cell selected!", type = "error")
            output$omim_link <- renderUI(NULL)  # Hide the link
          }
        })
}

shinyApp(ui = ui, server = server)

