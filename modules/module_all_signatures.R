# module_all_signatures.R

# UI do módulo
mod_all_signatures_ui <- function(id, dataset) {
  ns <- NS(id)  # Namespace to evitar conflitos de ID
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h4("Select Nomenclature"),
        selectizeInput(ns("signatureFilter"), "1. Select Signature:", choices = NULL, multiple = FALSE),
        hr(),
        br(), br(),
        downloadButton(ns("downloadData"), "Download Summary Results", icon = icon("download"))
      ),
      mainPanel(
        h3("Integrative Summary"),
        uiOutput(ns("summaryText")),  
        h3("Table Results"),
        tableOutput(ns("selectedRowTable"))  
      )
    )
  )
}


# Server do módulo
mod_all_signatures_server <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns  # Namespace
    
    # **Atualizar as opções do selectizeInput com todas as Signatures disponíveis**
    observe({
      req(dataset())
      data <- dataset()
      
      updateSelectizeInput(session, "signatureFilter", 
                           choices = unique(na.omit(data$Nomenclature)), 
                           selected = NULL, 
                           server = TRUE)  # 🔹 Enable server-side processing
    })
    
    # **Dados filtrados apenas pela Nomenclature selecionada**
    filteredData <- reactive({
      req(dataset(), input$signatureFilter)
      data <- dataset()
      filtered <- data[data$Nomenclature == input$signatureFilter, ]
      
      if (nrow(filtered) > 0) {
        return(filtered[1, , drop = FALSE])  
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
      filename = function() { paste0("selected_summary_", Sys.Date(), ".csv") },
      content = function(file) {
        data <- filteredData()
        data <- data %>% 
          mutate(Summary = paste0("Signature ", data$Nomenclature, " summary..."))
        write.csv(data, file, row.names = FALSE)
      }
    )
    
  })
}

