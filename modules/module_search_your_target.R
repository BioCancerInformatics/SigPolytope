# # UI Module
# mod_search_your_target_ui <- function(id) {
#   ns <- NS(id)  # Namespace to avoid ID conflicts
# 
#   tagList(
#     h3("Search Your Gene"),
#     textInput(ns("targetFilter"), "Enter Gene:", placeholder = "Type gene name here..."),
#     uiOutput(ns("targetStatus"), style = "margin-top: 10px;"),
#     DTOutput(ns("resultsTable")),
#     downloadButton(ns("downloadData"), "Download Results", style = "margin-top: 20px;")
#   )
# }
# 
# # Server Module
# mod_search_your_target_server <- function(id, dataset) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
# 
#     # Check if the manually entered gene exists in the Gene column
#     target_exists <- reactive({
#       req(dataset(), input$targetFilter != "")
#       data <- dataset()
#       gene_input <- trimws(toupper(input$targetFilter))  # Standardize input
#       target_found <- any(toupper(data$Target) == gene_input)
#       return(target_found)
#     })
# 
#     # Display the status of the target gene
#     output$targetStatus <- renderUI({
#       if (input$targetFilter == "") {
#         return(NULL)
#       }
# 
#       if (!target_exists()) {
#         return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
#       } else {
#         return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
#       }
#     })
#     
#     # Filter Target based on the target gene in the Signature column
#     filtered_data <- reactive({
#       req(input$targetFilter != "", target_exists())
#       target_pattern <- paste0("\\b", trimws(toupper(input$targetFilter)), "\\b")
#       filtered <- Target[grepl(target_pattern, toupper(Target$Signature)), ]
#       return(filtered)
#     })
# 
#     # Render the results table with the filtered data
#     output$resultsTable <- renderDT({
#       req(filtered_data())
# 
#       data <- filtered_data()
# 
#       # Torna a coluna Nomenclature clicável e realmente copiável
#       data$Nomenclature <- sprintf(
#         '<span style="cursor: pointer; color: #007bff;"
#             onclick="navigator.clipboard.writeText(`%s`).then(() => {
#             console.log(`Copied: %s`);
#             }).catch(err => console.error(`Copy failed`, err));">%s</span>',
#         data$Nomenclature, data$Nomenclature, data$Nomenclature
#       )
# 
#       datatable(data,
#                 escape = FALSE,  # Permite renderização de HTML na coluna Nomenclature
#                 options = list(pageLength = 10, autoWidth = TRUE),
#                 rownames = FALSE)
#     })
# 
#     # Enable download of the filtered data
#     output$downloadData <- downloadHandler(
#       filename = function() {
#         paste0("Filtered_Results_", input$targetFilter, ".csv")
#       },
#       content = function(file) {
#         data_to_download <- filtered_data()
#         write.csv(data_to_download, file, row.names = FALSE)
#       }
#     )
# 
#   })
# }


# mod_search_your_target_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     h3("Search Your Gene"),
#     textInput(ns("targetFilter"), "Enter Gene:", placeholder = "Type gene name here..."),
#     uiOutput(ns("targetStatus"), style = "margin-top: 10px;"),
#     DTOutput(ns("resultsTable")),
#     downloadButton(ns("downloadData"), "Download Results", style = "margin-top: 20px;"),
#     
#     tags$hr(),  # Separator
#     
#     uiOutput(ns("signatureList"), style = "margin-top: 20px;")  # Signature list
#   )
# }
# 
# mod_search_your_target_server <- function(id, search_dataset, Target_dataset) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     # Check if gene exists in search_dataset
#     target_exists <- reactive({
#       req(search_dataset(), input$targetFilter != "")
#       data <- search_dataset()
#       gene_input <- trimws(toupper(input$targetFilter))
#       any(toupper(data$Target) == gene_input)
#     })
#     
#     # Status output
#     output$targetStatus <- renderUI({
#       if (input$targetFilter == "") {
#         return(NULL)
#       }
#       if (!target_exists()) {
#         return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
#       } else {
#         return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
#       }
#     })
#     
#     # Filtered data for resultsTable
#     filtered_data <- reactive({
#       req(input$targetFilter != "", target_exists())
#       gene_input <- trimws(toupper(input$targetFilter))
#       
#       data <- search_dataset()
#       pattern <- paste0("\\b", gene_input, "\\b")
#       
#       filtered <- data[grepl(pattern, toupper(data$Signature)), ]
#       return(filtered)
#     })
#     
#     # Render filtered results table
#     output$resultsTable <- renderDT({
#       req(filtered_data())
#       
#       data <- filtered_data()
#       
#       data$Nomenclature <- sprintf(
#         '<span style="cursor: pointer; color: #007bff;"
#             onclick="navigator.clipboard.writeText(`%s`).then(() => {
#             console.log(`Copied: %s`);
#             }).catch(err => console.error(`Copy failed`, err));">%s</span>',
#         data$Nomenclature, data$Nomenclature, data$Nomenclature
#       )
#       
#       datatable(data,
#                 escape = FALSE,
#                 options = list(pageLength = 10, autoWidth = TRUE),
#                 rownames = FALSE)
#     })
#     
#     # Download handler
#     output$downloadData <- downloadHandler(
#       filename = function() {
#         paste0("Filtered_Results_", input$targetFilter, ".csv")
#       },
#       content = function(file) {
#         data_to_download <- filtered_data()
#         write.csv(data_to_download, file, row.names = FALSE)
#       }
#     )
#     
#     # Render signature list + corresponding nomenclature from Target_dataset
#     output$signatureList <- renderUI({
#       req(input$targetFilter != "", Target_dataset())
#       gene_input <- trimws(toupper(input$targetFilter))
#       
#       # Corrected: Target_dataset()$Signatures
#       match_vector <- sapply(Target_dataset()$Signatures, function(sig) {
#         sig_clean <- gsub("[()]", "", sig)
#         genes <- unlist(strsplit(sig_clean, "\\+"))
#         genes <- trimws(genes)
#         any(toupper(genes) == gene_input)
#       })
#       
#       filtered_target <- Target_dataset()[match_vector, ]
#       
#       if (nrow(filtered_target) == 0) {
#         return(HTML("<p><em>No signatures contain this gene.</em></p>"))
#       }
#       
#       signatures <- filtered_target$Signatures
#       nomenclatures <- filtered_target$Nomenclature
#       
#       sig_items <- paste0("- ", signatures, collapse = "<br>")
#       nom_items <- paste0("- ", nomenclatures, collapse = "<br>")
#       
#       tags$div(
#         h4("Signatures containing this gene:"),
#         HTML(sig_items),
#         br(), br(),
#         h4("Corresponding Nomenclature:"),
#         HTML(nom_items)
#       )
#     })
#   })
# }

# mod_search_your_target_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     h3("Search Your Gene"),
#     
#     fluidRow(
#       column(6,
#              textInput(ns("targetFilter"), "Enter Gene:", placeholder = "Type gene name here..."),
#              uiOutput(ns("targetStatus"), style = "margin-top: 10px;"),
#              shinycssloaders::withSpinner(DT::DTOutput(ns("resultsTable")), type = 6),
#              downloadButton(ns("downloadData"), "Download Results", style = "margin-top: 20px;")
#       ),
#       
#       column(6,
#              wellPanel(
#                tags$h4("Corresponding Nomenclature:"),
#                div(
#                  shinycssloaders::withSpinner(uiOutput(ns("nomenclatureList")), type = 6),
#                  style = "max-height: 300px; overflow-y: auto; border: 1px solid #ccc; border-radius: 8px; padding: 10px;"
#                ),
#                tags$hr(),
#                tags$h4("Corresponding Signature:"),
#                shinycssloaders::withSpinner(uiOutput(ns("correspondingSignature")), type = 6),
#                tags$hr(),
#                tags$h4("Nomenclature Interpretation:"),
#                shinycssloaders::withSpinner(uiOutput(ns("nomenclatureInterpretation")), type = 6),
#                downloadButton(ns("downloadInterpretation"), "Download Interpretation", style = "margin-top: 10px;")
#              )
#       )
#     )
#   )
# }




# mod_search_your_target_server <- function(id, search_dataset, Target_dataset) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
# 
#     # Store selected Nomenclature
#     selected_nom <- reactiveVal(NULL)
# 
#     # Check if gene exists
#     target_exists <- reactive({
#       req(search_dataset(), input$targetFilter != "")
#       data <- search_dataset()
#       gene_input <- trimws(toupper(input$targetFilter))
#       any(toupper(data$Target) == gene_input)
#     })
# 
#     # Immediate Status output
#     output$targetStatus <- renderUI({
#       if (input$targetFilter == "") return(NULL)
#       if (!target_exists()) {
#         return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
#       } else {
#         return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
#       }
#     })
# 
#     # Filtered data for resultsTable
#     filtered_data <- reactive({
#       req(input$targetFilter != "", target_exists())
#       gene_input <- trimws(toupper(input$targetFilter))
#       data <- search_dataset()
#       pattern <- paste0("\\b", gene_input, "\\b")
#       data[grepl(pattern, toupper(data$Signature)), ]
#     })
# 
#     output$resultsTable <- renderDT({
#       req(filtered_data())
#       data <- filtered_data()
#       data$Nomenclature <- sprintf(
#         '<span style="cursor: pointer; color: #007bff;"
#           onclick="navigator.clipboard.writeText(`%s`).then(() => {
#           console.log(`Copied: %s`);
#           }).catch(err => console.error(`Copy failed`, err));">%s</span>',
#         data$Nomenclature, data$Nomenclature, data$Nomenclature
#       )
#       datatable(data, escape = FALSE, options = list(pageLength = 10, autoWidth = TRUE), rownames = FALSE)
#     })
# 
#     output$downloadData <- downloadHandler(
#       filename = function() paste0("Filtered_Results_", input$targetFilter, ".csv"),
#       content = function(file) {
#         write.csv(filtered_data(), file, row.names = FALSE)
#       }
#     )
# 
#     # Event: when input changes, lookup once
#     matched_data <- eventReactive(input$targetFilter, {
#       req(input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
# 
#       match_vector <- sapply(Target_dataset()$Signatures, function(sig) {
#         sig_clean <- gsub("[()]", "", sig)
#         genes <- unlist(strsplit(sig_clean, "\\+"))
#         genes <- trimws(genes)
#         any(toupper(genes) == gene_input)
#       })
# 
#       Target_dataset()[match_vector, ]
#     }, ignoreNULL = TRUE, ignoreInit = TRUE)
# 
#     # Render Nomenclature clickable list
#     output$nomenclatureList <- renderUI({
#       data <- matched_data()
#       if (nrow(data) == 0) {
#         return(HTML("<p><em>No signatures contain this gene.</em></p>"))
#       }
# 
#       tags$div(
#         h4("Corresponding Nomenclature:"),
#         lapply(seq_len(nrow(data)), function(i) {
#           actionLink(ns(paste0("nom_", i)), label = data$Nomenclature[i], style = "display: block;")
#         })
#       )
#     })
# 
#     # Observe clicks on all generated actionLinks
#     observe({
#       data <- matched_data()
#       lapply(seq_len(nrow(data)), function(i) {
#         observeEvent(input[[paste0("nom_", i)]], {
#           selected_nom(data$Nomenclature[i])
#         }, ignoreInit = TRUE)
#       })
#     })
# 
#     # Show corresponding Signature for selected Nomenclature
#     output$correspondingSignature <- renderUI({
#       req(selected_nom())
#       data <- matched_data()
#       idx <- which(data$Nomenclature == selected_nom())
#       if (length(idx) == 0) return(NULL)
#       tags$div(
#         h4("Corresponding Signature:"),
#         tags$p(data$Signatures[idx[1]])
#       )
#     })
#   })
# }

# mod_search_your_target_server <- function(id, search_dataset, Target_dataset, tcga_types) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     selected_nomenclature <- reactiveVal(NULL)
#     selected_signature <- reactiveVal(NULL)
#     
#     target_exists <- reactive({
#       req(search_dataset(), input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
#       any(toupper(search_dataset()$Target) == gene_input)
#     })
#     
#     output$targetStatus <- renderUI({
#       if (input$targetFilter == "") return(NULL)
#       if (!target_exists()) {
#         HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>")
#       } else {
#         HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>")
#       }
#     })
#     
#     filtered_data <- reactive({
#       req(input$targetFilter != "", target_exists())
#       pattern <- paste0("\\b", trimws(toupper(input$targetFilter)), "\\b")
#       search_dataset()[grepl(pattern, toupper(search_dataset()$Signature)), ]
#     })
#     
#     output$resultsTable <- DT::renderDT({
#       req(filtered_data())
#       datatable(filtered_data(), options = list(pageLength = 10, autoWidth = TRUE), rownames = FALSE)
#     })
#     
#     output$downloadData <- downloadHandler(
#       filename = function() paste0("Filtered_Results_", input$targetFilter, ".csv"),
#       content = function(file) {
#         write.csv(filtered_data(), file, row.names = FALSE)
#       }
#     )
#     
#     output$nomenclatureList <- renderUI({
#       req(filtered_data())
#       nms <- filtered_data()$Nomenclature
#       tags$ul(
#         lapply(nms, function(nm) {
#           tags$li(
#             tags$a(nm, href = "#", onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("clicked_nomenclature"), nm))
#           )
#         })
#       )
#     })
#     
#     observeEvent(input$clicked_nomenclature, {
#       selected_nomenclature(input$clicked_nomenclature)
#       matched <- filtered_data()[filtered_data()$Nomenclature == input$clicked_nomenclature, ]
#       selected_signature(matched$Signature[1])
#     })
#     
#     output$correspondingSignature <- renderUI({
#       req(selected_signature())
#       HTML(paste0("<strong>", selected_signature(), "</strong>"))
#     })
#     
#     output$nomenclatureInterpretation <- renderUI({
#       req(selected_nomenclature())
#       interpret_nomenclature(selected_nomenclature(), tcga_types, Target_dataset())
#     })
#     
#     output$downloadInterpretation <- downloadHandler(
#       filename = function() paste0("Nomenclature_Interpretation_", selected_nomenclature(), ".txt"),
#       content = function(file) {
#         text <- gsub("<[^>]+>", "", interpret_nomenclature(selected_nomenclature(), tcga_types, Target_dataset()))
#         writeLines(text, file)
#       }
#     )
#     
#     interpret_nomenclature <- function(nomenclature, tcga_types, Target) {
#       components <- unlist(strsplit(nomenclature, "[.-]"))
#       if (length(components) < 11) return("Invalid nomenclature format.")
#       
#       cancer_abbr <- components[1]
#       cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
#       GSI <- components[2]
#       GFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
#       GFC <- GFC_codes[as.numeric(components[3])]
#       PFC_codes <- c("TMB", "MSI", "TSM")
#       PFC <- ifelse(components[4] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[4])], "NA")
#       SCS_codes <- c("Negative", "Positive")
#       SCS <- SCS_codes[ifelse(components[5] == "P", 2, 1)]
#       TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
#       TNC <- TNC_codes[as.numeric(components[6]) + 1]
#       
#       nomenclature_match <- Target[Target$Nomenclature == nomenclature, ]
#       
#       HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Type_Cox_", col)
#         if (!is.null(nomenclature_match[[col_name]])) paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
#         else paste0("<strong>", col, ":</strong> No data")
#       })
#       
#       SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Type_log_rank_", col)
#         if (!is.null(nomenclature_match[[col_name]])) paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
#         else paste0("<strong>", col, ":</strong> No data")
#       })
#       
#       TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
#       TMC <- TMC_codes[as.numeric(components[9])]
#       TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
#       TIC <- TIC_codes[as.numeric(components[10])]
#       RCD <- components[11]
#       
#       HTML(paste0(
#         "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
#         "<strong>GSI:</strong> ", GSI, "<br>",
#         "<strong>Multi-omic feature:</strong> ", GFC, "<br>",
#         "<strong>Phenotypic feature:</strong> ", PFC, "<br>",
#         "<strong>Spearman Sign:</strong> ", SCS, "<br>",
#         "<strong>Tumor Expression:</strong> ", TNC, "<br>",
#         "<strong>Hazard Ratio:</strong><br>", paste(HRC_values, collapse = "<br>"), "<br>",
#         "<strong>Kaplan-Meier:</strong><br>", paste(SMC_values, collapse = "<br>"), "<br>",
#         "<strong>TME Code:</strong> ", TMC, "<br>",
#         "<strong>TIL Code:</strong> ", TIC, "<br>",
#         "<strong>RCD forms:</strong> ", RCD
#       ))
#     }
#   })
# }

# mod_search_your_target_server <- function(id, search_dataset, Target_dataset) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     # Store selected Nomenclature
#     selected_nom <- reactiveVal(NULL)
#     
#     # Check if gene exists
#     target_exists <- reactive({
#       req(search_dataset(), input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
#       data <- search_dataset()
#       found <- any(toupper(data$Target) == gene_input)
#       cat("Gene input:", gene_input, "Exists:", found, "\n")
#       return(found)
#     })
#     
#     # Immediate Status output
#     output$targetStatus <- renderUI({
#       if (input$targetFilter == "") return(NULL)
#       if (!target_exists()) {
#         HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>")
#       } else {
#         HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>")
#       }
#     })
#     
#     # Filtered data for resultsTable
#     filtered_data <- reactive({
#       req(input$targetFilter != "", target_exists())
#       gene_input <- trimws(toupper(input$targetFilter))
#       data <- search_dataset()
#       
#       # Using robust composite signature parsing
#       pattern <- paste0("(^|[[:space:]+\\(])", gene_input, "($|[[:space:]+\\)])")
#       filtered <- data[grepl(pattern, toupper(as.character(data$Signature))), ]
#       
#       cat("Filtered Data Rows:", nrow(filtered), "\n")
#       return(filtered)
#     })
#     
#     output$resultsTable <- renderDT({
#       req(filtered_data())
#       data <- filtered_data()
#       
#       # Make Nomenclature clickable to copy
#       data$Nomenclature <- sprintf(
#         '<span style="cursor: pointer; color: #007bff;"
#           onclick="navigator.clipboard.writeText(`%s`).then(() => {
#           console.log(`Copied: %s`);
#           }).catch(err => console.error(`Copy failed`, err));">%s</span>',
#         data$Nomenclature, data$Nomenclature, data$Nomenclature
#       )
#       
#       datatable(data, escape = FALSE, options = list(pageLength = 10, autoWidth = TRUE), rownames = FALSE)
#     })
#     
#     output$downloadData <- downloadHandler(
#       filename = function() paste0("Filtered_Results_", input$targetFilter, ".csv"),
#       content = function(file) {
#         write.csv(filtered_data(), file, row.names = FALSE)
#       }
#     )
#     
#     # Event: when input changes, lookup once
#     matched_data <- eventReactive(input$targetFilter, {
#       req(input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
#       data <- Target_dataset()
#       
#       match_vector <- sapply(data$Signatures, function(sig) {
#         sig_clean <- gsub("[()]", "", sig)
#         genes <- unlist(strsplit(sig_clean, "\\+"))
#         genes <- trimws(genes)
#         any(toupper(genes) == gene_input)
#       })
#       
#       matched <- data[match_vector, ]
#       cat("Matched Nomenclatures:", nrow(matched), "\n")
#       return(matched)
#     }, ignoreNULL = TRUE, ignoreInit = TRUE)
#     
#     # Render Nomenclature clickable list
#     output$nomenclatureList <- renderUI({
#       data <- matched_data()
#       if (nrow(data) == 0) {
#         return(HTML("<p><em>No signatures contain this gene.</em></p>"))
#       }
#       
#       tags$div(
#         tags$h4("Corresponding Nomenclature:"),
#         lapply(seq_len(nrow(data)), function(i) {
#           actionLink(ns(paste0("nom_", i)), label = data$Nomenclature[i], style = "display: block;")
#         })
#       )
#     })
#     
#     # Observe clicks on all generated actionLinks
#     observe({
#       data <- matched_data()
#       lapply(seq_len(nrow(data)), function(i) {
#         observeEvent(input[[paste0("nom_", i)]], {
#           cat("Selected Nomenclature:", data$Nomenclature[i], "\n")
#           selected_nom(data$Nomenclature[i])
#         }, ignoreInit = TRUE)
#       })
#     })
#     
#     # Show corresponding Signature for selected Nomenclature
#     output$correspondingSignature <- renderUI({
#       req(selected_nom())
#       data <- matched_data()
#       idx <- which(data$Nomenclature == selected_nom())
#       if (length(idx) == 0) return(NULL)
#       tags$div(
#         tags$h4("Corresponding Signature:"),
#         tags$p(data$Signatures[idx[1]])
#       )
#     })
#   })
# }

# mod_search_your_target_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     h3("Search Your Gene"),
#     
#     fluidRow(
#       
#       # COLUNA ESQUERDA: Search + Signature
#       column(3,
#              wellPanel(
#                textInput(ns("targetFilter"), "Enter Gene:", placeholder = "Type gene name here..."),
#                actionButton(ns("searchGene"), "Search Gene"),
#                uiOutput(ns("targetStatus"), style = "margin-top: 10px;")
#              ),
#              br(),
#              wellPanel(
#                tags$h5("Corresponding Signature:"),
#                uiOutput(ns("correspondingSignature")),
#                br(),
#                downloadButton(ns("DownloadData"), "Download Data"),
#                style = "min-height: 200px; 
#                         border: 1px solid #ccc; 
#                         border-radius: 8px; 
#                         padding: 10px; 
#                         background-color: #f9f9f9; 
#                         box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       ),
#       
#       # COLUNA CENTRAL: Nomenclature
#       column(5,
#              wellPanel(
#                tags$h5("Corresponding Nomenclature:"),
#                tags$div(
#                  uiOutput(ns("nomenclatureList")),
#                  style = "max-height: 450px; 
#                     overflow-y: auto; 
#                     border: 1px solid #ccc; 
#                     border-radius: 8px; 
#                     padding: 10px; 
#                     background-color: #f9f9f9; 
#                     box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#                ),
#                style = "min-height: 450px; 
#                   border: 1px solid #ccc; 
#                   border-radius: 8px; 
#                   padding: 10px; 
#                   background-color: #f9f9f9; 
#                   box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       ),
#       
#       # COLUNA DIREITA: Interpretation
#       column(4,
#              wellPanel(
#                tags$h5("Interpretation:"),
#                shinycssloaders::withSpinner(uiOutput(ns("interpretation")), type = 4),
#                br(),
#                downloadButton(ns("downloadInterpretation"), "Download Interpretation"),
#                style = "min-height: 450px; 
#                         border: 1px solid #ccc; 
#                         border-radius: 8px; 
#                         padding: 10px; 
#                         background-color: #f9f9f9; 
#                         box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       )
#     )
#   )
# }
# 
# 
# 
# mod_search_your_target_server <- function(id, search_dataset, Target_dataset, tcga_types) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     selected_nom <- reactiveVal(NULL)
#     
#     target_exists <- eventReactive(input$searchGene, {
#       req(search_dataset(), input$targetFilter != "")
#       data <- search_dataset()
#       if (is.null(data) || nrow(data) == 0 || !"Target" %in% names(data)) return(FALSE)
#       gene_input <- trimws(toupper(input$targetFilter))
#       any(toupper(data$Target) == gene_input)
#     })
#     
#     output$targetStatus <- renderUI({
#       req(input$searchGene)
#       if (is.null(search_dataset()) || nrow(search_dataset()) == 0) {
#         return(HTML("<p><strong style='color: red;'>❌ Dataset not loaded.</strong></p>"))
#       }
#       if (!target_exists()) {
#         return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
#       } else {
#         return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
#       }
#     })
#     
#     matched_data <- eventReactive(input$searchGene, {
#       req(input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
#       target_data <- Target_dataset()
#       if (is.null(target_data) || nrow(target_data) == 0 || !"Signatures" %in% names(target_data)) return(data.frame())
#       signatures <- toupper(target_data$Signatures)
#       signatures[is.na(signatures)] <- ""
#       pattern <- paste0("(^|\\s|\\+|\\(|\\)|`)", gene_input, "($|\\s|\\+|\\(|\\)|`)")
#       match_vector <- grepl(pattern, signatures)
#       filtered <- target_data[match_vector, , drop = FALSE]
#       filtered
#     }, ignoreNULL = TRUE, ignoreInit = TRUE)
#     
#     output$nomenclatureList <- renderUI({
#       req(input$searchGene)
#       data <- matched_data()
#       
#       if (nrow(data) == 0) {
#         return(HTML("<p><em>No signatures contain this gene.</em></p>"))
#       }
#       
#       tagList(
#         # ----- CSS da lista -----
#         tags$style(HTML("
#       .nom-item {
#         display: flex;
#         align-items: center;
#         justify-content: space-between;
#         padding: 6px 8px;
#         margin-bottom: 6px;
#         border-radius: 6px;
#         background-color: #ffffff;
#         border: 1px solid #dcdcdc;
#       }
#       .nom-item:hover {
#         background-color: #eef3f8;
#       }
#       .nom-label {
#         font-family: 'Fira Code', 'Consolas', 'Courier New', monospace;
#         font-size: 12px;
#         word-break: break-all;
#         padding-right: 10px;
#         color: #333;
#         text-decoration: none;
#       }
#       .nom-label:hover {
#         text-decoration: underline;
#       }
#       .nom-copy-btn {
#         border: none;
#         background: transparent;
#         cursor: pointer;
#         color: #007bff;
#         padding: 0 6px;
#       }
#       .nom-copy-btn:hover {
#         color: #0056b3;
#       }
#     ")),
#         
#         # (o título principal já está no wellPanel, não repetimos aqui)
#         
#         tags$div(
#           lapply(seq_len(nrow(data)), function(i) {
#             nom <- data$Nomenclature[i]
#             
#             tags$div(
#               class = "nom-item",
#               
#               # texto da nomenclature – AGORA é um actionLink, como antes
#               actionLink(
#                 inputId = ns(paste0("nom_", i)),
#                 label   = nom,
#                 class   = "nom-label"
#               ),
#               
#               # botão copiar
#               tags$button(
#                 type  = "button",
#                 class = "nom-copy-btn",
#                 title = "Copy to clipboard",
#                 onclick = sprintf(
#                 "(function(txt){
#                 if (navigator.clipboard && navigator.clipboard.writeText) {
#                 navigator.clipboard.writeText(txt);
#      } else {
#      var ta = document.createElement('textarea');
#      ta.value = txt;
#      document.body.appendChild(ta);
#      ta.select();
#      document.execCommand('copy');
#      document.body.removeChild(ta);
#      }
#      })(%s);
#      this.blur();
#                 ",
#                   jsonlite::toJSON(nom, auto_unbox = TRUE)
#                 )
#                 ,
#                 shiny::icon("copy")
#               )
#             )
#           })
#         )
#       )
#     })
#     
#     observe({
#       data <- matched_data()
#       lapply(seq_len(nrow(data)), function(i) {
#         observeEvent(input[[paste0("nom_", i)]], {
#           selected_nom(data$Nomenclature[i])
#         }, ignoreInit = TRUE)
#       })
#     })
#     
#     output$correspondingSignature <- renderUI({
#       req(selected_nom())
#       data <- matched_data()
#       idx <- which(data$Nomenclature == selected_nom())
#       if (length(idx) == 0) return(NULL)
#       tags$div(
#         h4("Corresponding Signature:"),
#         tags$p(data$Signatures[idx[1]])
#       )
#     })
#     
#     output$interpretation <- renderUI({
#       req(selected_nom())
#       nomenclature <- selected_nom()
#       components <- unlist(strsplit(nomenclature, "[.-]"))
#       if (length(components) < 14) return(HTML("<p><em>Invalid nomenclature format.</em></p>"))
#       cancer_abbr <- components[1]
#       cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
#       USI <- components[2]
#       OFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
#       OFC <- OFC_codes[as.numeric(components[7])]
#       PFC_codes <- c("Tumor Mutation Burden", "Microsatellite instability", "Stemness")
#       PFC <- ifelse(components[8] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[8])], "NA")
#       SCS_codes <- c("Negative", "Positive")
#       SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
#       TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
#       TNC <- TNC_codes[as.numeric(components[10]) + 1]
#       
#       # Buscar a linha correspondente no data frame Target
#       nomenclature_match <- Target_dataset()[Target_dataset()$Nomenclature == nomenclature, ]
#       
#       # Extrair as novas informações
#       metabolism <- nomenclature_match$Metabolism
#       Pathways <- nomenclature_match$Pathways
#       interaction <- nomenclature_match$Meaningful_interaction
#       metabolic_cd <- nomenclature_match$Metabolic_cell_death
#       
#       HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Cox_", col, "_type")
#         if (!is.null(nomenclature_match[[col_name]])) paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]]) else paste0("<strong>", col, ":</strong> No data")
#       })
#       SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0(col, "_worst_prognosis_group")
#         if (!is.null(nomenclature_match[[col_name]])) paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]]) else paste0("<strong>", col, ":</strong> No data")
#       })
#       TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
#       TMC <- TMC_codes[as.numeric(components[13])]
#       TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
#       TIC <- TIC_codes[as.numeric(components[14])]
# 
#       interpretation <- paste0(
#         "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
#         "<strong>USI:</strong> ", USI, " - Unique signature identifier<br>",
#         "<strong>Metabolism:</strong> ", metabolism, "<br>",
#         "<strong>Pathway:</strong> ", Pathways, "<br>",
#         "<strong>Meaningful interaction:</strong> ", interaction, "<br>",
#         "<strong>Metabolic cell death:</strong> ", metabolic_cd, "<br>",
#         "<strong>Omic layer:</strong> ", components[7], " – ", OFC, "<br>",
#         "<strong>Phenotypic layer:</strong> ", components[8], " – ", PFC, "<br>",
#         "<strong>Spearman Correlation Sign:</strong> ", components[9], " – ", SCS, "<br>",
#         "<strong>Tumor vs. Non-Tumor Tissue:</strong> ", components[10], " – ", TNC, "<br>",
#         "<strong>Hazard Ratio:</strong> ", components[11], "<br>",
#         paste(HRC_values, collapse = "<br>"), "<br>",
#         "<strong>Kaplan-Meier survival classification:</strong>", components[12], "<br>",
#         paste(SMC_values, collapse = "<br>"), "<br>",
#         "<strong>Tumor Microenvironment:</strong> ", components[13], " – ", TMC, "<br>",
#         "<strong>Tumor Immune Infiltration:</strong> ", components[14], " – ", TIC, "<br>"
#       )
#       HTML(interpretation)
#     })
#     
# 
#     
#     generate_interpretation_text <- function(nomenclature, tcga_types, Target_dataset) {
#       components <- unlist(strsplit(nomenclature, "[.-]"))
#       if (length(components) < 14) return("Invalid nomenclature format.")
#       
#       cancer_abbr <- components[1]
#       cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
#       USI <- components[2]
#       
#       OFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
#       OFC <- OFC_codes[as.numeric(components[7])]
#       
#       PFC_codes <- c("TMB", "MSI", "TSM")
#       PFC <- ifelse(components[8] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[8])], "NA")
#       
#       SCS_codes <- c("Negative", "Positive")
#       SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
#       
#       TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
#       TNC <- TNC_codes[as.numeric(components[10]) + 1]
#       
#       nomenclature_match <- Target_dataset()[Target_dataset()$Nomenclature == nomenclature, ]
#       
#       # Extrair as novas informações
#       metabolism <- nomenclature_match$Metabolism
#       Pathways <- nomenclature_match$Pathways
#       interaction <- nomenclature_match$Meaningful_interaction
#       metabolic_cd <- nomenclature_match$Metabolic_cell_death
#       
#       HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Cox_", col, "_type")
#         if (!is.null(nomenclature_match[[col_name]])) paste0(col, ": ", nomenclature_match[[col_name]]) else paste0(col, ": No data")
#       })
#       
#       SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0(col, "_worst_prognosis_group")
#         if (!is.null(nomenclature_match[[col_name]])) paste0(col, ": ", nomenclature_match[[col_name]]) else paste0(col, ": No data")
#       })
#       
#       TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
#       TMC <- TMC_codes[as.numeric(components[13])]
#       
#       TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
#       TIC <- TIC_codes[as.numeric(components[14])]
#       
#       interpretation <- paste0(
#         "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
#         "<strong>USI:</strong> ", USI, " - Unique signature identifier<br>",
#         "<strong>Metabolism:</strong> ", metabolism, "<br>",
#         "<strong>Pathway:</strong> ", Pathways, "<br>",
#         "<strong>Meaningful interaction:</strong> ", interaction, "<br>",
#         "<strong>Metabolic cell death:</strong> ", metabolic_cd, "<br>",
#         "<strong>Omic layer:</strong> ", components[7], " – ", OFC, "<br>",
#         "<strong>Phenotypic layer:</strong> ", components[8], " – ", PFC, "<br>",
#         "<strong>Spearman Correlation Sign:</strong> ", components[9], " – ", SCS, "<br>",
#         "<strong>Tumor vs. Non-Tumor Tissue:</strong> ", components[10], " – ", TNC, "<br>",
#         "<strong>Hazard Ratio:</strong> ", components[11], "<br>",
#         paste(HRC_values, collapse = "<br>"), "<br>",
#         "<strong>Kaplan-Meier survival classification:</strong>", components[11], "<br>",
#         paste(SMC_values, collapse = "<br>"), "<br>",
#         "<strong>Tumor Microenvironment:</strong> ", components[13], " – ", TMC, "<br>",
#         "<strong>Tumor Immune Infiltration:</strong> ", components[14], " – ", TIC, "<br>"
#       )
#       
#       return(interpretation)
#     }
#     
#     
#     output$downloadInterpretation <- downloadHandler(
#       filename = function() paste0("Interpretation_", selected_nom(), ".txt"),
#       content = function(file) {
#         req(selected_nom())
#         interpretation <- generate_interpretation_text(selected_nom(), tcga_types, Target_dataset())
#         writeLines(interpretation, con = file)
#       }
#     )
#     
#     output$downloadDataRow <- downloadHandler(
#       filename = function() paste0("Data_Row_", selected_nom(), ".csv"),
#       content = function(file) {
#         req(selected_nom())
#         data <- Target_dataset()
#         row <- data[data$Nomenclature == selected_nom(), , drop = FALSE]
#         write.csv(row, file, row.names = FALSE)
#       }
#     )
#     
#   })
# }



# mod_search_your_target_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     h3("Search Your Gene"),
#     
#     fluidRow(
#       
#       # COLUNA ESQUERDA: Search + Signature
#       column(3,
#              wellPanel(
#                textInput(ns("targetFilter"), "Enter Gene:", placeholder = "Type gene name here..."),
#                actionButton(ns("searchGene"), "Search Gene"),
#                uiOutput(ns("targetStatus"), style = "margin-top: 10px;")
#              ),
#              br(),
#              wellPanel(
#                tags$h5("Corresponding Signature:"),
#                uiOutput(ns("correspondingSignature")),
#                br(),
#                # 🔁 corrigido para bater com o output$downloadDataRow no server
#                downloadButton(ns("downloadDataRow"), "Download Data"),
#                style = "min-height: 200px; 
#                         border: 1px solid #ccc; 
#                         border-radius: 8px; 
#                         padding: 10px; 
#                         background-color: #f9f9f9; 
#                         box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       ),
#       
#       # COLUNA CENTRAL: Nomenclature
#       column(5,
#              wellPanel(
#                tags$h5("Corresponding Nomenclature:"),
#                tags$div(
#                  uiOutput(ns("nomenclatureList")),
#                  style = "max-height: 450px; 
#                     overflow-y: auto; 
#                     border: 1px solid #ccc; 
#                     border-radius: 8px; 
#                     padding: 10px; 
#                     background-color: #f9f9f9; 
#                     box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#                ),
#                style = "min-height: 450px; 
#                   border: 1px solid #ccc; 
#                   border-radius: 8px; 
#                   padding: 10px; 
#                   background-color: #f9f9f9; 
#                   box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       ),
#       
#       # COLUNA DIREITA: Interpretation + Regulatory network
#       column(4,
#              wellPanel(
#                tags$h5("Interpretation:"),
#                shinycssloaders::withSpinner(uiOutput(ns("interpretation")), type = 4),
#                br(),
#                downloadButton(ns("downloadInterpretation"), "Download Interpretation"),
#                hr(),
#                tags$h5("Regulatory circuitry network"),
#                # Mensagem se não houver rede
#                uiOutput(ns("regulatoryStatus")),
#                # Plot com spinner; só aparece se houver rede
#                shinycssloaders::withSpinner(
#                  plotOutput(ns("regulatoryPlot"), height = "350px"),
#                  type = 4
#                ),
#                style = "min-height: 450px; 
#                         border: 1px solid #ccc; 
#                         border-radius: 8px; 
#                         padding: 10px; 
#                         background-color: #f9f9f9; 
#                         box-shadow: 0 2px 5px rgba(0,0,0,0.1);"
#              )
#       )
#     )
#   )
# }
# 
# 
# mod_search_your_target_server <- function(id,
#                                           search_dataset,
#                                           Target_dataset,
#                                           tcga_types,
#                                           regulatory_dataset) {
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     selected_nom <- reactiveVal(NULL)
#     
#     # ----------------------------------------------------------
#     # 1) Checagem se o gene existe no search_dataset
#     # ----------------------------------------------------------
#     target_exists <- eventReactive(input$searchGene, {
#       req(search_dataset(), input$targetFilter != "")
#       data <- search_dataset()
#       if (is.null(data) || nrow(data) == 0 || !"Target" %in% names(data)) return(FALSE)
#       gene_input <- trimws(toupper(input$targetFilter))
#       any(toupper(data$Target) == gene_input)
#     })
#     
#     output$targetStatus <- renderUI({
#       req(input$searchGene)
#       if (is.null(search_dataset()) || nrow(search_dataset()) == 0) {
#         return(HTML("<p><strong style='color: red;'>❌ Dataset not loaded.</strong></p>"))
#       }
#       if (!target_exists()) {
#         return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
#       } else {
#         return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
#       }
#     })
#     
#     # ----------------------------------------------------------
#     # 2) Filtrar as signatures/nomenclatures do Target_dataset
#     # ----------------------------------------------------------
#     matched_data <- eventReactive(input$searchGene, {
#       req(input$targetFilter != "")
#       gene_input <- trimws(toupper(input$targetFilter))
#       target_data <- Target_dataset()
#       if (is.null(target_data) || nrow(target_data) == 0 || !"Signatures" %in% names(target_data)) {
#         return(data.frame())
#       }
#       
#       signatures <- toupper(target_data$Signatures)
#       signatures[is.na(signatures)] <- ""
#       pattern <- paste0("(^|\\s|\\+|\\(|\\)|`)", gene_input, "($|\\s|\\+|\\(|\\)|`)")
#       match_vector <- grepl(pattern, signatures)
#       filtered <- target_data[match_vector, , drop = FALSE]
#       filtered
#     }, ignoreNULL = TRUE, ignoreInit = TRUE)
#     
#     # ----------------------------------------------------------
#     # 3) Lista de nomenclatures clicáveis
#     # ----------------------------------------------------------
#     output$nomenclatureList <- renderUI({
#       req(input$searchGene)
#       data <- matched_data()
#       
#       if (nrow(data) == 0) {
#         return(HTML("<p><em>No signatures contain this gene.</em></p>"))
#       }
#       
#       tagList(
#         # CSS
#         tags$style(HTML("
#           .nom-item {
#             display: flex;
#             align-items: center;
#             justify-content: space-between;
#             padding: 6px 8px;
#             margin-bottom: 6px;
#             border-radius: 6px;
#             background-color: #ffffff;
#             border: 1px solid #dcdcdc;
#           }
#           .nom-item:hover {
#             background-color: #eef3f8;
#           }
#           .nom-label {
#             font-family: 'Fira Code', 'Consolas', 'Courier New', monospace;
#             font-size: 12px;
#             word-break: break-all;
#             padding-right: 10px;
#             color: #333;
#             text-decoration: none;
#           }
#           .nom-label:hover {
#             text-decoration: underline;
#           }
#           .nom-copy-btn {
#             border: none;
#             background: transparent;
#             cursor: pointer;
#             color: #007bff;
#             padding: 0 6px;
#           }
#           .nom-copy-btn:hover {
#             color: #0056b3;
#           }
#         ")),
#         
#         tags$div(
#           lapply(seq_len(nrow(data)), function(i) {
#             nom <- data$Nomenclature[i]
#             
#             tags$div(
#               class = "nom-item",
#               
#               # clique na nomenclature
#               actionLink(
#                 inputId = ns(paste0("nom_", i)),
#                 label   = nom,
#                 class   = "nom-label"
#               ),
#               
#               # botão copiar
#               tags$button(
#                 type  = "button",
#                 class = "nom-copy-btn",
#                 title = "Copy to clipboard",
#                 onclick = sprintf(
#                   "(function(txt){
#                     if (navigator.clipboard && navigator.clipboard.writeText) {
#                       navigator.clipboard.writeText(txt);
#                     } else {
#                       var ta = document.createElement('textarea');
#                       ta.value = txt;
#                       document.body.appendChild(ta);
#                       ta.select();
#                       document.execCommand('copy');
#                       document.body.removeChild(ta);
#                     }
#                   })(%s);
#                   this.blur();",
#                   jsonlite::toJSON(nom, auto_unbox = TRUE)
#                 ),
#                 shiny::icon("copy")
#               )
#             )
#           })
#         )
#       )
#     })
#     
#     # Observers para cada nomenclature clicada
#     observe({
#       data <- matched_data()
#       lapply(seq_len(nrow(data)), function(i) {
#         observeEvent(input[[paste0("nom_", i)]], {
#           selected_nom(data$Nomenclature[i])
#         }, ignoreInit = TRUE)
#       })
#     })
#     
#     # ----------------------------------------------------------
#     # 4) Corresponding Signature (do Target_dataset)
#     # ----------------------------------------------------------
#     output$correspondingSignature <- renderUI({
#       req(selected_nom())
#       data <- matched_data()
#       idx <- which(data$Nomenclature == selected_nom())
#       if (length(idx) == 0) return(NULL)
#       tags$div(
#         h4("Corresponding Signature:"),
#         tags$p(data$Signatures[idx[1]])
#       )
#     })
#     
#     # ----------------------------------------------------------
#     # 5) Interpretation textual (igual ao que você já tinha)
#     # ----------------------------------------------------------
#     output$interpretation <- renderUI({
#       req(selected_nom())
#       nomenclature <- selected_nom()
#       components <- unlist(strsplit(nomenclature, "[.-]"))
#       if (length(components) < 14) {
#         return(HTML("<p><em>Invalid nomenclature format.</em></p>"))
#       }
#       
#       cancer_abbr <- components[1]
#       cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
#       USI <- components[2]
#       
#       OFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
#       OFC <- OFC_codes[as.numeric(components[7])]
#       
#       PFC_codes <- c("Tumor Mutation Burden", "Microsatellite instability", "Stemness")
#       PFC <- ifelse(components[8] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[8])], "NA")
#       
#       SCS_codes <- c("Negative", "Positive")
#       SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
#       
#       TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
#       TNC <- TNC_codes[as.numeric(components[10]) + 1]
#       
#       nomenclature_match <- Target_dataset()[Target_dataset()$Nomenclature == nomenclature, ]
#       
#       metabolism    <- nomenclature_match$Metabolism
#       Pathways      <- nomenclature_match$Pathways
#       interaction   <- nomenclature_match$Meaningful_interaction
#       metabolic_cd  <- nomenclature_match$Metabolic_cell_death
#       
#       HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Cox_", col, "_type")
#         if (!is.null(nomenclature_match[[col_name]])) {
#           paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
#         } else {
#           paste0("<strong>", col, ":</strong> No data")
#         }
#       })
#       
#       SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0(col, "_worst_prognosis_group")
#         if (!is.null(nomenclature_match[[col_name]])) {
#           paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
#         } else {
#           paste0("<strong>", col, ":</strong> No data")
#         }
#       })
#       
#       TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
#       TMC <- TMC_codes[as.numeric(components[13])]
#       
#       TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
#       TIC <- TIC_codes[as.numeric(components[14])]
#       
#       interpretation <- paste0(
#         "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
#         "<strong>USI:</strong> ", USI, " - Unique signature identifier<br>",
#         "<strong>Metabolism:</strong> ", metabolism, "<br>",
#         "<strong>Pathway:</strong> ", Pathways, "<br>",
#         "<strong>Meaningful interaction:</strong> ", interaction, "<br>",
#         "<strong>Metabolic cell death:</strong> ", metabolic_cd, "<br>",
#         "<strong>Omic layer:</strong> ", components[7], " – ", OFC, "<br>",
#         "<strong>Phenotypic layer:</strong> ", components[8], " – ", PFC, "<br>",
#         "<strong>Spearman Correlation Sign:</strong> ", components[9], " – ", SCS, "<br>",
#         "<strong>Tumor vs. Non-Tumor Tissue:</strong> ", components[10], " – ", TNC, "<br>",
#         "<strong>Hazard Ratio:</strong> ", components[11], "<br>",
#         paste(HRC_values, collapse = "<br>"), "<br>",
#         "<strong>Kaplan-Meier survival classification:</strong>", components[12], "<br>",
#         paste(SMC_values, collapse = "<br>"), "<br>",
#         "<strong>Tumor Microenvironment:</strong> ", components[13], " – ", TMC, "<br>",
#         "<strong>Tumor Immune Infiltration:</strong> ", components[14], " – ", TIC, "<br>"
#       )
#       
#       HTML(interpretation)
#     })
#     
#     # ----------------------------------------------------------
#     # 6) Função para exportar a interpretação em texto
#     # ----------------------------------------------------------
#     generate_interpretation_text <- function(nomenclature, tcga_types, Target_dataset) {
#       components <- unlist(strsplit(nomenclature, "[.-]"))
#       if (length(components) < 14) return("Invalid nomenclature format.")
#       
#       cancer_abbr <- components[1]
#       cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
#       USI <- components[2]
#       
#       OFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
#       OFC <- OFC_codes[as.numeric(components[7])]
#       
#       PFC_codes <- c("TMB", "MSI", "TSM")
#       PFC <- ifelse(components[8] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[8])], "NA")
#       
#       SCS_codes <- c("Negative", "Positive")
#       SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
#       
#       TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
#       TNC <- TNC_codes[as.numeric(components[10]) + 1]
#       
#       nomenclature_match <- Target_dataset[Target_dataset$Nomenclature == nomenclature, , drop = FALSE]
#       
#       metabolism    <- nomenclature_match$Metabolism
#       Pathways      <- nomenclature_match$Pathways
#       interaction   <- nomenclature_match$Meaningful_interaction
#       metabolic_cd  <- nomenclature_match$Metabolic_cell_death
#       
#       HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0("Cox_", col, "_type")
#         if (!is.null(nomenclature_match[[col_name]])) {
#           paste0(col, ": ", nomenclature_match[[col_name]])
#         } else {
#           paste0(col, ": No data")
#         }
#       })
#       
#       SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
#         col_name <- paste0(col, "_worst_prognosis_group")
#         if (!is.null(nomenclature_match[[col_name]])) {
#           paste0(col, ": ", nomenclature_match[[col_name]])
#         } else {
#           paste0(col, ": No data")
#         }
#       })
#       
#       TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
#       TMC <- TMC_codes[as.numeric(components[13])]
#       
#       TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
#       TIC <- TIC_codes[as.numeric(components[14])]
#       
#       interpretation <- paste0(
#         cancer_abbr, " – ", cancer_full, "\n",
#         "USI: ", USI, " - Unique signature identifier\n",
#         "Metabolism: ", metabolism, "\n",
#         "Pathway: ", Pathways, "\n",
#         "Meaningful interaction: ", interaction, "\n",
#         "Metabolic cell death: ", metabolic_cd, "\n",
#         "Omic layer: ", components[7], " – ", OFC, "\n",
#         "Phenotypic layer: ", components[8], " – ", PFC, "\n",
#         "Spearman Correlation Sign: ", components[9], " – ", SCS, "\n",
#         "Tumor vs. Non-Tumor Tissue: ", components[10], " – ", TNC, "\n",
#         "Hazard Ratio: ", components[11], "\n",
#         paste(HRC_values, collapse = "\n"), "\n",
#         "Kaplan-Meier survival classification: ", components[11], "\n",
#         paste(SMC_values, collapse = "\n"), "\n",
#         "Tumor Microenvironment: ", components[13], " – ", TMC, "\n",
#         "Tumor Immune Infiltration: ", components[14], " – ", TIC, "\n"
#       )
#       
#       interpretation
#     }
#     
#     output$downloadInterpretation <- downloadHandler(
#       filename = function() paste0("Interpretation_", selected_nom(), ".txt"),
#       content = function(file) {
#         req(selected_nom())
#         interpretation <- generate_interpretation_text(selected_nom(), tcga_types, Target_dataset())
#         writeLines(interpretation, con = file)
#       }
#     )
#     
#     # ----------------------------------------------------------
#     # 7) Download da linha de dados (Target_dataset)
#     # ----------------------------------------------------------
#     output$downloadDataRow <- downloadHandler(
#       filename = function() paste0("Data_Row_", selected_nom(), ".csv"),
#       content = function(file) {
#         req(selected_nom())
#         data <- Target_dataset()
#         row <- data[data$Nomenclature == selected_nom(), , drop = FALSE]
#         write.csv(row, file, row.names = FALSE)
#       }
#     )
#     
#     # ----------------------------------------------------------
#     # 8) Integração com regulatory_dataset + grafo
#     # ----------------------------------------------------------
#     
#     # Filtra regulatory_dataset pela nomenclature selecionada
#     reg_matches <- reactive({
#       req(selected_nom())
#       df_reg <- regulatory_dataset()
#       
#       if (is.null(df_reg) ||
#           nrow(df_reg) == 0 ||
#           !"Nomenclature" %in% names(df_reg)) {
#         return(NULL)
#       }
#       
#       hits <- df_reg[df_reg$Nomenclature == selected_nom(), , drop = FALSE]
#       if (nrow(hits) == 0) return(NULL)
#       
#       hits
#     })
#     
#     # Mensagem sobre existência ou não de rede
#     output$regulatoryStatus <- renderUI({
#       req(selected_nom())
#       hits <- reg_matches()
#       if (is.null(hits)) {
#         HTML("<p><em>There is no regulatory circuitry available for this nomenclature.</em></p>")
#       } else {
#         HTML("<p><strong>Regulatory circuitry available for this nomenclature.</strong></p>")
#       }
#     })
#     
#     # Plot da rede (só quando há matches)
#     output$regulatoryPlot <- renderPlot({
#       hits <- reg_matches()
#       validate(
#         need(!is.null(hits) && nrow(hits) > 0,
#              "")  # mensagem vazia porque o texto já vai no regulatoryStatus
#       )
#       net_obj <- build_network_from_matches(hits)
#       net_obj$plot
#     }, res = 120)
#     
#   })
# }
# 



mod_search_your_target_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    
    # ---------- ESTILO GLOBAL ----------
    tags$head(
      tags$style(HTML("
        .syg-hero {
          background: #f4f6f9;
          padding: 20px;
          border-radius: 10px;
          margin-bottom: 20px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }
        .syg-card {
          background: white;
          border-radius: 12px;
          padding: 18px;
          box-shadow: 0 2px 6px rgba(0,0,0,0.06);
          border: 1px solid #e8e8e8;
          margin-bottom: 20px;
        }
        /* mesma altura nas três colunas */
        .syg-card-fixed {
          min-height: 520px;
          display: flex;
          flex-direction: column;
        }
        .syg-title {
          font-size: 18px;
          font-weight: 600;
          margin-bottom: 12px;
        }
        .syg-nom-list {
          flex: 1;
          overflow-y: auto;
          padding-right: 4px;
          margin-top: 8px;
        }
        .syg-nom-pill {
          display: block;
          width: 100%;
          text-align: left;
          padding: 6px 10px;
          margin-bottom: 4px;
          background: #f5f7fb;
          border-radius: 6px;
          border: 1px solid #dde2ee;
          font-size: 11px;
          font-family: 'Fira Code','Consolas','Courier New',monospace;
          cursor: pointer;
          color: #333;
          white-space: nowrap;
          text-overflow: ellipsis;
          overflow: hidden;
        }
        .syg-nom-pill:hover {
          background: #e3e7f3;
        }
        .syg-nom-pager {
          margin-top: 8px;
          font-size: 11px;
          color: #666;
        }
      "))
    ),
    
    # ---------- BARRA DE BUSCA ----------
    div(
      class = "syg-hero",
      h3("Search Your Target"),
      p("Type a gene, miRNA or lncRNA to explore signatures and regulatory networks."),
      textInput(ns("targetFilter"), label = NULL,
                placeholder = "e.g. FADS1, hsa-miR-224-5p, LINC00511"),
      actionButton(ns("searchGene"), "Search Target"),
      br(), br(),
      uiOutput(ns("targetStatus"))
    ),
    
    # ---------- LAYOUT DOS RESULTADOS ----------
    fluidRow(
      
      # COLUNA 1 – Nomenclatures
      column(
        width = 3,
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("searchGene")),
          div(
            class = "syg-card syg-card-fixed",
            div(class = "syg-title", "Nomenclatures"),
            div(class = "syg-nom-list",
                uiOutput(ns("nomenclatureList"))
            ),
            div(class = "syg-nom-pager",
                uiOutput(ns("nomenclaturePager"))
            )
          )
        )
      ),
      
      # COLUNA 2 – Interpretation
      column(
        width = 3,
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("searchGene")),
          div(
            class = "syg-card syg-card-fixed",
            div(class = "syg-title", "Interpretation"),
            div(
              style = "flex:1; overflow-y:auto; padding-right:4px;",
              shinycssloaders::withSpinner(
                uiOutput(ns("interpretation")),
                type = 4
              )
            ),
            hr(),
            downloadButton(ns("downloadDataRow"),
                           "Download Data",
                           style = "width:100%; margin-bottom:8px;"),
            downloadButton(ns("downloadInterpretation"),
                           "Download Interpretation",
                           style = "width:100%;")
          )
        )
      ),
      
      # COLUNA 3 – Regulatory circuitry
      column(
        width = 6,
        conditionalPanel(
          condition = sprintf("input['%s'] > 0", ns("searchGene")),
          div(
            class = "syg-card syg-card-fixed",
            div(class = "syg-title", "Regulatory networks"),
            uiOutput(ns("regulatoryStatus")),
            div(
              style = "flex:1;",
              shinycssloaders::withSpinner(
                plotOutput(ns("regulatoryPlot"), height = "500px"),
                type = 4
              )
            ),
            hr(),
            downloadButton(ns("downloadPlot"),
                           "Download Plot (PDF)",
                           style = "width:100%;")
          )
        )
      )
    )
  )
}







mod_search_your_target_server <- function(id,
                                          search_dataset,
                                          Target_dataset,
                                          tcga_types,
                                          regulatory_dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    selected_nom <- reactiveVal(NULL)
    
    # --- controle de página do "carrossel" de nomenclaturas ---
    page_nom  <- reactiveVal(1L)
    page_size <- 15L  # nº de nomenclaturas por página
    
    # ----------------------------------------------------------
    # 1) Checagem se o gene existe no search_dataset
    # ----------------------------------------------------------
    target_exists <- eventReactive(input$searchGene, {
      req(search_dataset(), input$targetFilter != "")
      data <- search_dataset()
      if (is.null(data) || nrow(data) == 0 || !"Target" %in% names(data)) return(FALSE)
      gene_input <- trimws(toupper(input$targetFilter))
      any(toupper(data$Target) == gene_input)
    })
    
    output$targetStatus <- renderUI({
      req(input$searchGene)
      if (is.null(search_dataset()) || nrow(search_dataset()) == 0) {
        return(HTML("<p><strong style='color: red;'>❌ Dataset not loaded.</strong></p>"))
      }
      if (!target_exists()) {
        return(HTML("<p><strong style='color: red;'>❌ This gene is not in our database.</strong></p>"))
      } else {
        return(HTML("<p><strong style='color: green;'>✅ This gene is in our database. Go ahead!</strong></p>"))
      }
    })
    
    # ----------------------------------------------------------
    # 2) Filtrar as signatures/nomenclatures do Target_dataset
    # ----------------------------------------------------------
    matched_data <- eventReactive(input$searchGene, {
      req(input$targetFilter != "")
      gene_input <- trimws(toupper(input$targetFilter))
      target_data <- Target_dataset()
      if (is.null(target_data) || nrow(target_data) == 0 || !"Signatures" %in% names(target_data)) {
        return(data.frame())
      }
      
      signatures <- toupper(target_data$Signatures)
      signatures[is.na(signatures)] <- ""
      pattern <- paste0("(^|\\s|\\+|\\(|\\)|`)", gene_input, "($|\\s|\\+|\\(|\\)|`)")
      match_vector <- grepl(pattern, signatures)
      filtered <- target_data[match_vector, , drop = FALSE]
      filtered
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    # sempre que buscar um novo alvo, volta para a página 1
    observeEvent(input$searchGene, {
      page_nom(1L)
    })
    
    # ----------------------------------------------------------
    # 3) Lista de nomenclatures clicáveis (com paginação)
    # ----------------------------------------------------------
    output$nomenclatureList <- renderUI({
      req(input$searchGene)
      data <- matched_data()
      
      if (is.null(data) || nrow(data) == 0) {
        return(HTML("<p><em>No signatures contain this target.</em></p>"))
      }
      
      n <- nrow(data)
      total_pages <- max(1L, ceiling(n / page_size))
      cur <- max(1L, min(page_nom(), total_pages))
      
      idx_start <- (cur - 1L) * page_size + 1L
      idx_end   <- min(cur * page_size, n)
      idx_seq   <- idx_start:idx_end
      
      tagList(
        lapply(idx_seq, function(i) {
          nom <- data$Nomenclature[i]
          # apenas o texto clicável, sem ícone de copiar
          actionLink(
            inputId = ns(paste0("nom_", i)),
            label   = nom,
            class   = "syg-nom-pill"
          )
        })
      )
    })
    
    # pager (Previous / Next + texto "Page X of Y")
    output$nomenclaturePager <- renderUI({
      data <- matched_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      n <- nrow(data)
      total_pages <- max(1L, ceiling(n / page_size))
      cur <- max(1L, min(page_nom(), total_pages))
      
      tagList(
        fluidRow(
          column(
            width = 4,
            actionButton(ns("nom_prev"), "Previous",
                         style = "width:100%; padding:2px 6px;")
          ),
          column(
            width = 4,
            div(style = "text-align:center; padding-top:4px;",
                sprintf("Page %d of %d", cur, total_pages))
          ),
          column(
            width = 4,
            actionButton(ns("nom_next"), "Next",
                         style = "width:100%; padding:2px 6px;")
          )
        )
      )
    })
    
    observeEvent(input$nom_prev, {
      page_nom(max(1L, page_nom() - 1L))
    })
    
    observeEvent(input$nom_next, {
      data <- matched_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      total_pages <- max(1L, ceiling(nrow(data) / page_size))
      page_nom(min(total_pages, page_nom() + 1L))
    })
    
    # Observers para cada nomenclature clicada
    observe({
      data <- matched_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      lapply(seq_len(nrow(data)), function(i) {
        observeEvent(input[[paste0("nom_", i)]], {
          selected_nom(data$Nomenclature[i])
        }, ignoreInit = TRUE)
      })
    })
    
    # ----------------------------------------------------------
    # 4) Interpretation textual (incluindo a Corresponding Signature)
    # ----------------------------------------------------------
    output$interpretation <- renderUI({
      req(selected_nom())
      nomenclature <- selected_nom()
      components <- unlist(strsplit(nomenclature, "[.-]"))
      if (length(components) < 14) {
        return(HTML("<p><em>Invalid nomenclature format.</em></p>"))
      }
      
      cancer_abbr <- components[1]
      cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
      USI <- components[2]
      
      OFC_codes <- c("Protein Expression", "Mutations", "CNV",
                     "miRNA Expression", "Transcript Expression",
                     "mRNA Expression", "CpG Methylation")
      OFC <- OFC_codes[as.numeric(components[7])]
      
      PFC_codes <- c("Tumor Mutation Burden", "Microsatellite instability", "Stemness")
      PFC <- ifelse(components[8] %in% c("1", "2", "3"),
                    PFC_codes[as.numeric(components[8])],
                    "NA")
      
      SCS_codes <- c("Negative", "Positive")
      SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
      
      TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
      TNC <- TNC_codes[as.numeric(components[10]) + 1]
      
      # Linha correspondente no Target_dataset
      target_all <- Target_dataset()
      nomenclature_match <- target_all[target_all$Nomenclature == nomenclature, , drop = FALSE]
      
      metabolism    <- nomenclature_match$Metabolism
      Pathways      <- nomenclature_match$Pathways
      interaction   <- nomenclature_match$Meaningful_interaction
      metabolic_cd  <- nomenclature_match$Metabolic_cell_death
      corr_signature <- nomenclature_match$Signatures[1]
      
      HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0("Cox_", col, "_type")
        if (!is.null(nomenclature_match[[col_name]])) {
          paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
        } else {
          paste0("<strong>", col, ":</strong> No data")
        }
      })
      
      SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0(col, "_worst_prognosis_group")
        if (!is.null(nomenclature_match[[col_name]])) {
          paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]])
        } else {
          paste0("<strong>", col, ":</strong> No data")
        }
      })
      
      TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
      TMC <- TMC_codes[as.numeric(components[13])]
      
      TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
      TIC <- TIC_codes[as.numeric(components[14])]
      
      interpretation <- paste0(
        "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
        "<strong>USI:</strong> ", USI, " - Unique signature identifier<br>",
        "<strong>Corresponding signature:</strong> ", corr_signature, "<br><br>",
        "<strong>Metabolism:</strong> ", metabolism, "<br>",
        "<strong>Pathway:</strong> ", Pathways, "<br>",
        "<strong>Meaningful interaction:</strong> ", interaction, "<br>",
        "<strong>Metabolic cell death:</strong> ", metabolic_cd, "<br>",
        "<strong>Omic layer:</strong> ", components[7], " – ", OFC, "<br>",
        "<strong>Phenotypic layer:</strong> ", components[8], " – ", PFC, "<br>",
        "<strong>Spearman Correlation Sign:</strong> ", components[9], " – ", SCS, "<br>",
        "<strong>Tumor vs. Non-Tumor Tissue:</strong> ", components[10], " – ", TNC, "<br>",
        "<strong>Hazard Ratio:</strong> ", components[11], "<br>",
        paste(HRC_values, collapse = "<br>"), "<br>",
        "<strong>Kaplan-Meier survival classification:</strong> ", components[12], "<br>",
        paste(SMC_values, collapse = "<br>"), "<br>",
        "<strong>Tumor Microenvironment:</strong> ", components[13], " – ", TMC, "<br>",
        "<strong>Tumor Immune Infiltration:</strong> ", components[14], " – ", TIC, "<br>"
      )
      
      HTML(interpretation)
    })
    
    # ----------------------------------------------------------
    # 5) Função para exportar a interpretação em texto
    # ----------------------------------------------------------
    generate_interpretation_text <- function(nomenclature, tcga_types, Target_dataset) {
      components <- unlist(strsplit(nomenclature, "[.-]"))
      if (length(components) < 14) return("Invalid nomenclature format.")
      
      cancer_abbr <- components[1]
      cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
      USI <- components[2]
      
      OFC_codes <- c("Protein Expression", "Mutations", "CNV",
                     "miRNA Expression", "Transcript Expression",
                     "mRNA Expression", "CpG Methylation")
      OFC <- OFC_codes[as.numeric(components[7])]
      
      PFC_codes <- c("TMB", "MSI", "TSM")
      PFC <- ifelse(components[8] %in% c("1", "2", "3"),
                    PFC_codes[as.numeric(components[8])],
                    "NA")
      
      SCS_codes <- c("Negative", "Positive")
      SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
      
      TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
      TNC <- TNC_codes[as.numeric(components[10]) + 1]
      
      nomenclature_match <- Target_dataset[Target_dataset$Nomenclature == nomenclature, , drop = FALSE]
      
      metabolism    <- nomenclature_match$Metabolism
      Pathways      <- nomenclature_match$Pathways
      interaction   <- nomenclature_match$Meaningful_interaction
      metabolic_cd  <- nomenclature_match$Metabolic_cell_death
      corr_signature <- nomenclature_match$Signatures[1]
      
      HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0("Cox_", col, "_type")
        if (!is.null(nomenclature_match[[col_name]])) {
          paste0(col, ": ", nomenclature_match[[col_name]])
        } else {
          paste0(col, ": No data")
        }
      })
      
      SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0(col, "_worst_prognosis_group")
        if (!is.null(nomenclature_match[[col_name]])) {
          paste0(col, ": ", nomenclature_match[[col_name]])
        } else {
          paste0(col, ": No data")
        }
      })
      
      TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
      TMC <- TMC_codes[as.numeric(components[13])]
      
      TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
      TIC <- TIC_codes[as.numeric(components[14])]
      
      paste0(
        cancer_abbr, " – ", cancer_full, "\n",
        "USI: ", USI, " - Unique signature identifier\n",
        "Corresponding signature: ", corr_signature, "\n\n",
        "Metabolism: ", metabolism, "\n",
        "Pathway: ", Pathways, "\n",
        "Regulatory network: ", interaction, "\n",
        "Metabolic cell death: ", metabolic_cd, "\n",
        "Omic layer: ", components[7], " – ", OFC, "\n",
        "Phenotypic layer: ", components[8], " – ", PFC, "\n",
        "Spearman Correlation Sign: ", components[9], " – ", SCS, "\n",
        "Tumor vs. Non-Tumor Tissue: ", components[10], " – ", TNC, "\n",
        "Hazard Ratio: ", components[11], "\n",
        paste(HRC_values, collapse = "\n"), "\n",
        "Kaplan-Meier survival classification: ", components[12], "\n",
        paste(SMC_values, collapse = "\n"), "\n",
        "Tumor Microenvironment: ", components[13], " – ", TMC, "\n",
        "Tumor Immune Infiltration: ", components[14], " – ", TIC, "\n"
      )
    }
    
    output$downloadInterpretation <- downloadHandler(
      filename = function() paste0("Interpretation_", selected_nom(), ".txt"),
      content = function(file) {
        req(selected_nom())
        interpretation <- generate_interpretation_text(selected_nom(), tcga_types, Target_dataset())
        writeLines(interpretation, con = file)
      }
    )
    
    # ----------------------------------------------------------
    # 6) Download da linha de dados (Target_dataset)
    # ----------------------------------------------------------
    output$downloadDataRow <- downloadHandler(
      filename = function() paste0("Data_Row_", selected_nom(), ".csv"),
      content = function(file) {
        req(selected_nom())
        data <- Target_dataset()
        row <- data[data$Nomenclature == selected_nom(), , drop = FALSE]
        write.csv(row, file, row.names = FALSE)
      }
    )
    
    # ================= DOWNLOAD DO PLOT EM PDF =====================
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("Regulatory_network_", selected_nom(), ".pdf")
      },
      content = function(file) {
        
        hits <- reg_matches()
        if (is.null(hits) || nrow(hits) == 0) {
          pdf(file)
          plot.new()
          text(0.5, 0.5, "No regulatory network available")
          dev.off()
          return()
        }
        
        net_obj <- build_network_from_matches(hits)
        
        pdf(file, width = 11, height = 8.5)   # tamanho ideal horizontal
        print(net_obj$plot)
        dev.off()
      }
    )
    
    # ----------------------------------------------------------
    # 7) Integração com regulatory_dataset + grafo
    # ----------------------------------------------------------
    reg_matches <- reactive({
      req(selected_nom())
      df_reg <- regulatory_dataset()
      
      if (is.null(df_reg) ||
          nrow(df_reg) == 0 ||
          !"Nomenclature" %in% names(df_reg)) {
        return(NULL)
      }
      
      hits <- df_reg[df_reg$Nomenclature == selected_nom(), , drop = FALSE]
      if (nrow(hits) == 0) return(NULL)
      
      hits
    })
    
    output$regulatoryStatus <- renderUI({
      req(selected_nom())
      hits <- reg_matches()
      if (is.null(hits)) {
        HTML("<p><em>There is no regulatory network available for this nomenclature.</em></p>")
      } else {
        HTML("<p><strong>Regulatory network available for this nomenclature.</strong></p>")
      }
    })
    
    output$regulatoryPlot <- renderPlot({
      hits <- reg_matches()
      validate(
        need(!is.null(hits) && nrow(hits) > 0, "")
      )
      net_obj <- build_network_from_matches(hits)
      net_obj$plot
    }, res = 120)
  })
}

