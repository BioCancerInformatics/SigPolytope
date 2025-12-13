# module_signature.R

# UI do módulo
mod_signature_ui <- function(id, signature_name) {
  ns <- NS(id)  # Namespace to prevent ID conflicts
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h4("Step-by-Step Filtering"),
        uiOutput(ns("omicLayerFilterUI")),  # Omic Layer é fixo, mas dinâmico
        selectInput(ns("cancerFilter"), "2. Select Cancer Type:", choices = NULL),
        selectInput(ns("MetabolismFilter"), "3. Select Metabolism:", choices = NULL),
        selectInput(ns("PathwayFilter"), "4. Select Pathway:", choices = NULL),
        selectInput(ns("MetabolicCellDeathFilter"), "5. Select Metabolic Cell Death:", choices = NULL),
        selectizeInput(ns("signatureFilter"), "6. Select Signature:", choices = NULL, multiple = FALSE),
        hr(),
        actionButton(ns("resetFilters"), "Reset Filters", icon = icon("redo")),
        br(), br(),
        downloadButton(ns("downloadData"), "Download Summary Results", icon = icon("download"))
      ),
      mainPanel(
        h3("Integrative Summary"),
        uiOutput(ns("summaryText")),  
        h3("Table results"),
        tableOutput(ns("selectedRowTable"))  
      )
    )
  )
}

# Server do módulo atualizado
mod_signature_server <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns  # Use namespace to avoid ID collisions in Shiny modules
    
    # Helper function to update selectInput choices with fallback message
    update_filter <- function(input_id, choices) {
      if (length(choices) == 0) {
        choices <- "No results available"
      }
      updateSelectInput(session, input_id, choices = choices, selected = NULL)
    }
    
    # Render fixed but dynamic Omic Layer input (based on dataset content)
    output$omicLayerFilterUI <- renderUI({
      req(dataset())
      data <- dataset()
      omic_layer_value <- unique(na.omit(data$Omic_layer))[1]
      selectInput(
        ns("omicLayerFilter"),
        "1. Select Omic Layer:",
        choices = omic_layer_value,
        selected = omic_layer_value
      )
    })
    
    # Initialize Cancer Types based on the full dataset
    observe({
      req(dataset())
      data <- dataset()
      update_filter("cancerFilter", unique(na.omit(data$CTAB)))
    })
    
    # Upon selecting Cancer Type, update Metabolism options and reset downstream filters
    observeEvent(input$cancerFilter, {
      req(dataset(), input$cancerFilter)
      data <- dataset()
      filtered <- data[data$CTAB == input$cancerFilter, ]
      
      update_filter("MetabolismFilter", unique(na.omit(filtered$Metabolism)))
      update_filter("PathwayFilter", NULL)
      update_filter("MetabolicCellDeathFilter", NULL)
      update_filter("signatureFilter", NULL)
    })
    
    # After selecting Metabolism, update Pathway options and reset downstream filters
    observeEvent(input$MetabolismFilter, {
      req(dataset(), input$cancerFilter, input$MetabolismFilter)
      data <- dataset()
      filtered <- data[
        data$CTAB == input$cancerFilter &
          data$Metabolism == input$MetabolismFilter, ]
      
      update_filter("PathwayFilter", unique(na.omit(filtered$Pathways)))
      update_filter("MetabolicCellDeathFilter", NULL)
      update_filter("signatureFilter", NULL)
    })
    
    # After selecting Pathways, update Metabolic Cell Death options
    observeEvent(input$PathwayFilter, {
      req(dataset(), input$cancerFilter, input$MetabolismFilter, input$PathwayFilter)
      data <- dataset()
      filtered <- data[
        data$CTAB == input$cancerFilter &
          data$Metabolism == input$MetabolismFilter &
          data$Pathways == input$PathwayFilter, ]
      
      update_filter("MetabolicCellDeathFilter", unique(na.omit(filtered$Metabolic_cell_death)))
      update_filter("signatureFilter", NULL)
    })
    
    # After selecting Metabolic Cell Death, update Signature options
    observeEvent(input$MetabolicCellDeathFilter, {
      req(dataset(), input$cancerFilter, input$MetabolismFilter,
          input$PathwayFilter, input$MetabolicCellDeathFilter)
      data <- dataset()
      filtered <- data[
        data$CTAB == input$cancerFilter &
          data$Metabolism == input$MetabolismFilter &
          data$Pathways == input$PathwayFilter &
          data$Metabolic_cell_death == input$MetabolicCellDeathFilter, ]
      
      updateSelectizeInput(session, "signatureFilter",
                           choices = unique(na.omit(filtered$Signatures)),
                           server = TRUE)  # Enable server-side loading
    })
    
    # Reactive object holding the filtered row (based on all filters selected)
    filteredData <- reactive({
      req(dataset(), input$cancerFilter, input$MetabolismFilter,
          input$PathwayFilter, input$MetabolicCellDeathFilter, input$signatureFilter)
      data <- dataset()
      filtered <- data[
        data$CTAB == input$cancerFilter &
          data$Metabolism == input$MetabolismFilter &
          data$Pathways == input$PathwayFilter &
          data$Metabolic_cell_death == input$MetabolicCellDeathFilter &
          data$Signatures == input$signatureFilter, ]
      
      if (nrow(filtered) > 0) {
        return(filtered[1, , drop = FALSE])  # Return only the first matching row
      } else {
        return(NULL)
      }
    })
    
  
    
    
    
    
    
    # Display Summary without Bold Formatting
    output$summaryText <- renderUI({
      req(filteredData())  # Ensure data is available
      row <- filteredData()
      
      # Ensure row is not empty
      if (nrow(row) == 0) {
        return(HTML("<p><i>No data available for the selected signature.</i></p>"))
      }
      
      # Helper function to format statistical results safely, ensuring "NS" is removed
      format_stat <- function(type, pval) {
        if (!is.null(type) && !is.na(type) && type != "NS") { 
          return(paste0(type, " (p = ", sprintf("%.3g", as.numeric(pval)), ")"))
        }
        return(NULL)  # Exclude NS values
      }
      
      # Lookup cancer full name from tcga_types
      cancer_abbr <- row[["CTAB"]]  # This is the abbreviation (e.g., "UVM")
      cancer_full <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
      
      # If no match found, keep the abbreviation
      if (is.na(cancer_full) || cancer_full == "") {
        cancer_full <- cancer_abbr
      }
      
    
      
      # Tumor_vs_normal and correlation text logic
      rho <- row[["Correlation_rho"]]
      p_value <- row[["Correlation_p.adj"]]
      phenotype <- row[["Phenotypic_layer"]]
      omic_layer <- row[["Omic_layer"]]
      tumor_vs_normal <- trimws(row[["Tumor_vs_normal"]])  # remove espaços
      
      formatted_rho <- if (!is.na(rho)) sprintf("%.3f", as.numeric(rho)) else NA
      formatted_p <- if (!is.na(p_value)) formatC(as.numeric(p_value), format = "e", digits = 2) else NA
      
      # Frase com ou sem expressão tumoral
      if (!is.na(tumor_vs_normal) && tumor_vs_normal != "No data") {
        correlation_text <- paste0(
          "<p>This signature exhibits a correlation between <b>", omic_layer, "</b> and <b>", phenotype,
          "</b> (ρ = ", formatted_rho, ", p.adj = ", formatted_p, 
          "), and its expression remains <b>", tumor_vs_normal, 
          "</b> in tumor tissue compared to non-tumor tissue.</p>"
        )
      } else {
        correlation_text <- paste0(
          "<p>This signature exhibits a correlation between <b>", omic_layer, "</b> and <b>", phenotype,
          "</b> (ρ = ", formatted_rho, ", p.adj = ", formatted_p, ").</p>"
        )
      }
      
      
      
      
      # Cox Regression Summary (Filter out NS values)
      cox_results <- Filter(Negate(is.null), c(
        if (!is.null(row[["Cox_OS_type"]]) && row[["Cox_OS_type"]] != "NS") 
          paste0("OS: ", format_stat(row[["Cox_OS_type"]], row[["Cox_OS_p.value"]])),
        if (!is.null(row[["Cox_DSS_type"]]) && row[["Cox_DSS_type"]] != "NS") 
          paste0("DSS: ", format_stat(row[["Cox_DSS_type"]], row[["Cox_DSS_p.value"]])),
        if (!is.null(row[["Cox_DFI_type"]]) && row[["Cox_DFI_type"]] != "NS") 
          paste0("DFI: ", format_stat(row[["Cox_DFI_type"]], row[["Cox_DFI_p.value"]])),
        if (!is.null(row[["Cox_PFI_type"]]) && row[["Cox_PFI_type"]] != "NS") 
          paste0("PFI: ", format_stat(row[["Cox_PFI_type"]], row[["Cox_PFI_p.value"]]))
      ))
      cox_summary <- if (length(cox_results) > 0) paste(cox_results, collapse = ", ") 
      else "not significant in protection or risk."
      
      
      # Survival Analysis Summary (Filter out NS values)
      survival_results <- Filter(Negate(is.null), c(
        if (!is.null(row[["OS_worst_prognosis_group"]]) && row[["OS_worst_prognosis_group"]] != "NS") 
          paste0("OS: ", format_stat(row[["OS_worst_prognosis_group"]], row[["OS_p.value"]])),
        if (!is.null(row[["DSS_worst_prognosis_group"]]) && row[["DSS_worst_prognosis_group"]] != "NS") 
          paste0("DSS: ", format_stat(row[["DSS_worst_prognosis_group"]], row[["DSS_p.value"]])),
        if (!is.null(row[["DFI_worst_prognosis_group"]]) && row[["DFI_worst_prognosis_group"]] != "NS") 
          paste0("DFI: ", format_stat(row[["DFI_worst_prognosis_group"]], row[["DFI_p.value"]])),
        if (!is.null(row[["PFI_worst_prognosis_group"]]) && row[["PFI_worst_prognosis_group"]] != "NS") 
          paste0("PFI: ", format_stat(row[["PFI_worst_prognosis_group"]], row[["PFI_p.value"]]))
      ))
      survival_summary <- if (length(survival_results) > 0) paste(survival_results, collapse = ", ")
      else "not significant results in prognostic metrics."
      
      
      # Microenvironment and Immune Classification (robust logic)
      microenv <- trimws(row[["Microenvironment_classification"]])
      immune <- trimws(row[["Immune_classification"]])
      
      microenv_valid <- !is.na(microenv) && microenv != "NS"
      immune_valid <- !is.na(immune) && immune != "NS"
      
      if (microenv_valid && immune_valid) {
        microimmune_text <- paste0(
          "This signature is associated with a <b>", microenv,
          "</b> profile and an immune phenotype classified as <b>", immune, "</b>."
        )
      } else if (microenv_valid) {
        microimmune_text <- paste0(
          "This signature is associated with a <b>", microenv, "</b> profile."
        )
      } else if (immune_valid) {
        microimmune_text <- paste0(
          "This signature exhibits an immune phenotype classified as <b>", immune, "</b>."
        )
      } else {
        microimmune_text <- ""
      }
      
      # Metabolic Cell Death and Interaction logic (evita duplicação da interação)
      mcd <- row[["Metabolic_cell_death"]]
      interaction <- trimws(row[["Meaningful_interaction"]])
      mcd_interaction_text <- ""
      
      # Só gerar texto se MCD for relevante
      if (!is.na(mcd) && mcd != "Unrelated") {
        mcd_interaction_text <- paste0(
          "This signature is associated with <b>", mcd, "</b>."
        )
      }
      
      # Text about metabolism, pathway, and molecular class
      class_text <- paste0(
        "<p>The nomenclature <b>", row[["Nomenclature"]], "</b> represents a molecular signature identified in <b>",
        cancer_full, "</b>. Its component(s) belong to the <b>", row[["Molecular_class"]], "</b> category and are involved in <b>",
        row[["Metabolism"]], "</b> and <b>", row[["Pathways"]], "</b>.</p>"
      )
      
      # Significância molecular interpretável (nova frase logo após a primeira)
      molecular_class <- row[["Molecular_class"]]
      interaction <- trimws(row[["Meaningful_interaction"]])
      
      if (!is.na(interaction) && interaction != "No meaningful interaction") {
        interaction_sentence <- paste0(
          "<p>This <b>", molecular_class,
          "</b> signature exhibits clinically meaningful molecular interactions with <b>", 
          interaction, "</b>.</p>"
        )
      } else {
        interaction_sentence <- ""
      }
      
      
      # Final HTML summary output
      HTML(paste0(
        class_text,
        interaction_sentence,
        correlation_text,
        "<p>", mcd_interaction_text, "</p>",
        "<p>Cox regression analysis suggests that the <b>", omic_layer, 
        "</b> of this signature is associated with: <b>", cox_summary, "</b>.</p>",
        "<p>Survival analysis shows <b>", survival_summary, "</b>.</p>",
        "<p>", microimmune_text, "</p>",
        "<p><i>For detailed statistical results, click 'Download'. To visualize the results, copy the signature nomenclature, open the 'Analysis and Plotting' dropdown menu, and select the desired analysis.</i></p>"
      ))
    })
    
    
  
    
    # Display filtered data in table
    output$selectedRowTable <- renderTable({
      req(filteredData())
      filteredData()
    })
    
    # Download handler
    output$downloadData <- downloadHandler(
      filename = function() { paste0("selected_summary_", Sys.Date(), ".tsv") },
      content = function(file) {
        data <- filteredData()
        data <- data %>% 
          mutate(Summary = HTML(paste0("<p>Signature ", data$Nomenclature, " summary...</p>")))
        readr::write_tsv(data, file, row.names = FALSE)
      }
    )
    
    
    # **Reset filters properly**
    observeEvent(input$resetFilters, {
      updateSelectInput(session, "cancerFilter", choices = unique(na.omit(dataset()$CTAB)), selected = NULL)
      updateSelectInput(session, "MetabolismFilter", choices = NULL, selected = NULL)
      updateSelectInput(session, "PathwayFilter", choices = NULL, selected = NULL)
      updateSelectInput(session, "MetabolicCellDeathFilter", choices = NULL, selected = NULL)
      updateSelectInput(session, "signatureFilter", choices = NULL, selected = NULL)
    })
  })
}
