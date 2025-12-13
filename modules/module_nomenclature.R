# UI Module
mod_nomenclature_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Nomenclature Interpreter"),
    textInput(ns("nomenclatureInput"), "Enter Nomenclature:", placeholder = "Type nomenclature here..."),
    uiOutput(ns("nomenclatureStatus"), style = "margin-top: 10px;"),
    uiOutput(ns("interpretation")),  # Changed from verbatimTextOutput to uiOutput
    downloadButton(ns("downloadData"), "Download Interpretation", style = "margin-top: 20px;")
  )
}


# Server Module
mod_nomenclature_server <- function(id, tcga_types, Target) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    nomenclature_exists <- reactive({
      req(Target(), input$nomenclatureInput != "")
      nomenclature_input <- trimws(toupper(input$nomenclatureInput))
      any(toupper(Target()$Nomenclature) == nomenclature_input)
    })
    
    output$nomenclatureStatus <- renderUI({
      if (input$nomenclatureInput == "") {
        return(NULL)
      }
      if (!nomenclature_exists()) {
        return(HTML("<p><strong style='color: red;'>❌ This nomenclature is not in our database.</strong></p>"))
      } else {
        return(HTML("<p><strong style='color: green;'>✅ Nomenclature found. Interpretation available.</strong></p>"))
      }
    })
    
    interpret_nomenclature <- function(nomenclature, tcga_types, Target) {
      components <- unlist(strsplit(nomenclature, "[.-]"))
      if (length(components) < 11) return("Invalid nomenclature format.")
      
      cancer_abbr <- components[1]
      cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
      USI <- components[2]
      
      OFC_codes <- c("Protein Expression", "Mutations", "CNV", "miRNA Expression", "Transcript Expression", "mRNA Expression", "CpG Methylation")
      OFC <- OFC_codes[as.numeric(components[7])]
      
      PFC_codes <- c("Tumor Mutation Burden", "Microsatellite instability", "Stemness")
      PFC <- ifelse(components[8] %in% c("1", "2", "3"), PFC_codes[as.numeric(components[8])], "NA")
      
      SCS_codes <- c("Negative", "Positive")
      SCS <- SCS_codes[ifelse(components[9] == "P", 2, 1)]
      
      TNC_codes <- c("No data", "Unchanged", "Underexpressed", "Overexpressed")
      TNC <- TNC_codes[as.numeric(components[10]) + 1]
      
      # Buscar a linha correspondente no data frame Target
      nomenclature_match <- Target()[Target()$Nomenclature == nomenclature, ]
      
      # Garantir que a linha foi encontrada
      if (nrow(nomenclature_match) == 0) return("Nomenclature not found in Target.")
      
      # Extrair as novas informações
      metabolism <- nomenclature_match$Metabolism
      Pathways <- nomenclature_match$Pathways
      interaction <- nomenclature_match$Meaningful_interaction
      metabolic_cd <- nomenclature_match$Metabolic_cell_death
      
      # Valores de Cox
      HRC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0("Cox_", col, "_type")
        if (!is.null(nomenclature_match[[col_name]]) && !is.na(nomenclature_match[[col_name]])) {
          return(paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]]))
        } else {
          return(paste0("<strong>", col, ":</strong> No data"))
        }
      })
      
      
      # Valores de Kaplan-Meier
      SMC_values <- sapply(c("DSS", "DFI", "PFI", "OS"), function(col) {
        col_name <- paste0(col, "_worst_prognosis_group")
        if (!is.null(nomenclature_match[[col_name]])) {
          return(paste0("<strong>", col, ":</strong> ", nomenclature_match[[col_name]]))
        } else {
          return(paste0("<strong>", col, ":</strong> No data"))
        }
      })
      
      TMC_codes <- c("Anti-tumoral", "Dual", "Pro-tumoral", "No significant data")
      TMC <- TMC_codes[as.numeric(components[13])]
      
      TIC_codes <- c("Hot", "Variable", "Cold", "No significant data")
      TIC <- TIC_codes[as.numeric(components[14])]
      
      interpretation <- paste0(
        "<strong>", cancer_abbr, "</strong> – ", cancer_full, "<br>",
        "<strong>USI:</strong> ", USI, " - Unique signature identifier<br>",
        "<strong>Metabolism:</strong> ", metabolism, "<br>",
        "<strong>Pathway:</strong> ", Pathways, "<br>",
        "<strong>Meaningful interaction:</strong> ", interaction, "<br>",
        "<strong>Metabolic cell death:</strong> ", metabolic_cd, "<br>",
        "<strong>Omic layer:</strong> ", components[7], " – ", OFC, "<br>",
        "<strong>Phenotypic layer:</strong> ", components[8], " – ", PFC, "<br>",
        "<strong>Spearman Correlation Sign:</strong> ", components[9], " – ", SCS, "<br>",
        "<strong>Tumor vs. Non-Tumor Tissue:</strong> ", components[10], " – ", TNC, "<br>",
        "<strong>Hazard Ratio:</strong><br>", paste(HRC_values, collapse = "<br>"), "<br>",
        "<strong>Kaplan-Meier survival classification:</strong><br>", paste(SMC_values, collapse = "<br>"), "<br>",
        "<strong>Tumor Microenvironment:</strong> ", components[13], " – ", TMC, "<br>",
        "<strong>Tumor Immune Infiltration:</strong> ", components[14], " – ", TIC, "<br>"
      )
      
      return(HTML(interpretation))  # Return HTML instead of plain text
    }
    
    output$interpretation <- renderUI({
      req(nomenclature_exists())
      interpret_nomenclature(input$nomenclatureInput, tcga_types, Target())
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("Nomenclature_Interpretation_", input$nomenclatureInput, ".txt")
      },
      content = function(file) {
        writeLines(gsub("<[^>]+>", "", interpret_nomenclature(input$nomenclatureInput, tcga_types, Target())), file)
      }
    )
  })
}
