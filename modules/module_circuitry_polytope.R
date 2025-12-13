# # =====================================================================
# # modules/module_circuitry_polytope.R
# # Módulo: filtros → seleção de Circuitries_id → plot 3D + resumo
# # =====================================================================
# 
# suppressPackageStartupMessages({
#   library(shiny)
#   library(plotly)
#   library(dplyr)
#   library(stringr)
#   library(htmltools)
#   library(htmlwidgets)
# })
# 
# # ---------------------------------------------------------------------
# # UI do módulo
# # ---------------------------------------------------------------------
# mod_circuitry_polytope_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     h2("Geometric Regulatory Circuitries – 3D Explorer"),
#     
#     fluidRow(
#       # ===================== COLUNA ESQUERDA: FILTROS ==================
#       box(
#         width = 3,
#         title = "Filters (paper-driven geometry)",
#         status = "primary",
#         solidHeader = TRUE,
#         
#         # 1) CTAB – tipo de câncer
#         selectInput(
#           ns("fil_ctab"),
#           "Cancer type (CTAB):",
#           choices = NULL,        # preenchido no server
#           multiple = TRUE
#         ),
#         
#         # 2) Metabolic superfamily
#         selectInput(
#           ns("fil_metab"),
#           "Metabolic superfamily (Metabolism):",
#           choices = NULL,        # preenchido no server
#           multiple = TRUE
#         ),
#         
#         # 3) Metabolic pathway
#         selectInput(
#           ns("fil_path"),
#           "Metabolic pathway (Pathways):",
#           choices = NULL,        # preenchido no server
#           multiple = TRUE
#         ),
#         
#         # 4) Metabolic cell death
#         selectInput(
#           ns("fil_mcd"),
#           "Metabolic cell death:",
#           choices = NULL,        # preenchido no server
#           multiple = TRUE
#         ),
#         
#         tags$hr(),
#         # 5) Regime de discordância (barycenter_distance)
#         selectInput(
#           ns("fil_dist_imp"),
#           "Barycenter-distance regime:",
#           choices = NULL,        # distance_implication
#           multiple = TRUE
#         ),
#         
#         # 6) Regime de volume / complexidade (vol_implication)
#         selectInput(
#           ns("fil_vol_imp"),
#           "Volume / complexity regime:",
#           choices = NULL,        # vol_implication
#           multiple = TRUE
#         ),
#         
#         # 7) Faixa quantitativa de barycenter_distance
#         sliderInput(
#           ns("fil_bary"),
#           "Barycenter distance range:",
#           min   = 0,
#           max   = 1,
#           value = c(0, 1),
#           step  = 0.01
#         ),
#         
#         # 8) Faixa quantitativa de vol_ratio
#         sliderInput(
#           ns("fil_volratio"),
#           "Volume ratio (sig/int) range:",
#           min   = 0,
#           max   = 1,
#           value = c(0, 1),
#           step  = 0.01
#         ),
#         
#         tags$hr(),
#         # 9) Seletor de Circuitries_id (dependente dos filtros)
#         uiOutput(ns("ui_circuitry_selector")),
#         
#         tags$hr(),
#         htmlOutput(ns("ui_count")),
#         
#         tags$hr(),
#         downloadButton(ns("dl_poly_html"), "Download 3D polytope (.html)")
#       ),
#       
#       # ===================== COLUNA DIREITA: PLOT + RESUMO ===============
#       box(
#         width = 9,
#         title = "3D dual polytope: signature × interaction",
#         status = "primary",
#         solidHeader = TRUE,
#         plotlyOutput(ns("poly_plot"), height = "600px"),
#         tags$hr(),
#         uiOutput(ns("summary_ui"))
#       )
#     )
#   )
# }
# 
# # ---------------------------------------------------------------------
# # SERVER do módulo
# # tensor: objeto com tensor$meta (data.frame com colunas listadas)
# # embedding: matriz/data.frame de coordenadas latentes (usado por
# #            build_circuitry_polytope() e plot_circuitry_polytope()).
# # ---------------------------------------------------------------------
# mod_circuitry_polytope_server <- function(id, tensor, embedding) {
#   
#   moduleServer(id, function(input, output, session) {
#     ns <- session$ns
#     
#     # ---------------- META DE CIRCUITRIES ------------------------------
#     meta <- reactive({
#       validate(need(!is.null(tensor), "Tensor is NULL"))
#       validate(need(!is.null(tensor$meta), "tensor$meta missing"))
#       tensor$meta
#     })
#     
#     # ---------------- INICIALIZAR FILTROS ------------------------------
#     observeEvent(meta(), {
#       m <- meta()
#       req(nrow(m) > 0)
#       
#       # CTAB
#       if ("CTAB" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_ctab",
#           choices = sort(unique(na.omit(as.character(m$CTAB))))
#         )
#       }
#       
#       # Metabolism
#       if ("Metabolism" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_metab",
#           choices = sort(unique(na.omit(as.character(m$Metabolism))))
#         )
#       }
#       
#       # Pathways
#       if ("Pathways" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_path",
#           choices = sort(unique(na.omit(as.character(m$Pathways))))
#         )
#       }
#       
#       # Metabolic_cell_death
#       if ("Metabolic_cell_death" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_mcd",
#           choices = sort(unique(na.omit(as.character(m$Metabolic_cell_death))))
#         )
#       }
#       
#       # distance_implication
#       if ("distance_implication" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_dist_imp",
#           choices = sort(unique(na.omit(as.character(m$distance_implication))))
#         )
#       }
#       
#       # vol_implication
#       if ("vol_implication" %in% names(m)) {
#         updateSelectInput(
#           session, "fil_vol_imp",
#           choices = sort(unique(na.omit(as.character(m$vol_implication))))
#         )
#       }
#       
#       # sliders dinâmicos para barycenter_distance e vol_ratio
#       if ("barycenter_distance" %in% names(m)) {
#         rng <- range(m$barycenter_distance, na.rm = TRUE)
#         if (is.finite(rng[1]) && is.finite(rng[2])) {
#           updateSliderInput(
#             session, "fil_bary",
#             min   = floor(rng[1] * 10) / 10,
#             max   = ceiling(rng[2] * 10) / 10,
#             value = rng
#           )
#         }
#       }
#       if ("vol_ratio" %in% names(m)) {
#         rng <- range(m$vol_ratio, na.rm = TRUE)
#         if (is.finite(rng[1]) && is.finite(rng[2])) {
#           updateSliderInput(
#             session, "fil_volratio",
#             min   = floor(rng[1] * 10) / 10,
#             max   = ceiling(rng[2] * 10) / 10,
#             value = rng
#           )
#         }
#       }
#     }, ignoreInit = FALSE)
#     
#     # ---------------- APLICAR FILTROS ------------------------------
#     filtered_meta <- reactive({
#       m <- meta()
#       req(nrow(m) > 0)
#       
#       # 1) CTAB
#       if (!is.null(input$fil_ctab) && length(input$fil_ctab) > 0) {
#         m <- m[m$CTAB %in% input$fil_ctab, , drop = FALSE]
#       }
#       
#       # 2) Metabolism
#       if (!is.null(input$fil_metab) && length(input$fil_metab) > 0) {
#         m <- m[m$Metabolism %in% input$fil_metab, , drop = FALSE]
#       }
#       
#       # 3) Pathways
#       if (!is.null(input$fil_path) && length(input$fil_path) > 0) {
#         m <- m[m$Pathways %in% input$fil_path, , drop = FALSE]
#       }
#       
#       # 4) Metabolic_cell_death
#       if (!is.null(input$fil_mcd) && length(input$fil_mcd) > 0) {
#         m <- m[m$Metabolic_cell_death %in% input$fil_mcd, , drop = FALSE]
#       }
#       
#       # 5) distance_implication
#       if (!is.null(input$fil_dist_imp) && length(input$fil_dist_imp) > 0 &&
#           "distance_implication" %in% names(m)) {
#         m <- m[m$distance_implication %in% input$fil_dist_imp, , drop = FALSE]
#       }
#       
#       # 6) vol_implication
#       if (!is.null(input$fil_vol_imp) && length(input$fil_vol_imp) > 0 &&
#           "vol_implication" %in% names(m)) {
#         m <- m[m$vol_implication %in% input$fil_vol_imp, , drop = FALSE]
#       }
#       
#       # 7) barycenter_distance (numérico)
#       if ("barycenter_distance" %in% names(m) && !is.null(input$fil_bary)) {
#         r <- input$fil_bary
#         m <- m[!is.na(m$barycenter_distance) &
#                  m$barycenter_distance >= r[1] &
#                  m$barycenter_distance <= r[2], , drop = FALSE]
#       }
#       
#       # 8) vol_ratio (numérico)
#       if ("vol_ratio" %in% names(m) && !is.null(input$fil_volratio)) {
#         r <- input$fil_volratio
#         m <- m[!is.na(m$vol_ratio) &
#                  m$vol_ratio >= r[1] &
#                  m$vol_ratio <= r[2], , drop = FALSE]
#       }
#       
#       m
#     })
#     
#     # ---------------- UI: seletor de Circuitries_id ---------------------
#     output$ui_circuitry_selector <- renderUI({
#       m <- filtered_meta()
#       
#       if (nrow(m) == 0) {
#         return(tags$p(
#           style = "color:red; margin-top:10px;",
#           "No circuitries match the selected filters."
#         ))
#       }
#       
#       selectizeInput(
#         ns("sel_circ"),
#         "Available Circuitries_id:",
#         choices  = m$Circuitries_id,
#         selected = m$Circuitries_id[1],
#         multiple = FALSE,
#         options  = list(placeholder = "Select one circuitry")
#       )
#     })
#     
#     # Contador de quantos circuitries sobram
#     output$ui_count <- renderUI({
#       m <- filtered_meta()
#       n <- nrow(m)
#       if (n == 0) {
#         HTML("<p><b>0</b> circuitries after filters.</p>")
#       } else {
#         HTML(glue::glue("<p><b>{n}</b> circuitries after filters.</p>"))
#       }
#     })
#     
#     # ---------------- Linha selecionada -------------------------------
#     selected_row <- reactive({
#       m <- filtered_meta()
#       req(nrow(m) > 0)
#       
#       cid <- input$sel_circ
#       if (is.null(cid) || !cid %in% m$Circuitries_id) {
#         return(m[1, , drop = FALSE])
#       }
#       m[m$Circuitries_id == cid, , drop = FALSE]
#     })
#     
#     # ---------------- Resumo textual (apenas colunas fornecidas) -------
#     output$summary_ui <- renderUI({
#       row <- selected_row()
#       req(nrow(row) == 1)
#       
#       pick <- function(nm) if (nm %in% names(row)) row[[nm]][1] else NA
#       
#       cid   <- pick("Circuitries_id")
#       ctab  <- pick("CTAB")
#       nom_s <- pick("Nomenclature_sig")
#       sig   <- pick("Signatures")
#       class_s <- pick("Molecular_class_sig")
#       
#       nom_i <- pick("Nomenclature_int")
#       inter <- pick("Interaction")
#       class_i <- pick("Molecular_class_int")
#       
#       metab <- pick("Metabolism")
#       path  <- pick("Pathways")
#       mcd   <- pick("Metabolic_cell_death")
#       
#       bary  <- pick("barycenter_distance")
#       dist_imp <- pick("distance_implication")
#       
#       sig_vol <- pick("sig_hull_vol")
#       int_vol <- pick("int_hull_vol")
#       vr      <- pick("vol_ratio")
#       vol_imp <- pick("vol_implication")
#       
#       HTML(paste0(
#         "<h4>Circuitry summary</h4>",
#         "<p><b>Circuitries_id:</b> ", cid, "<br>",
#         "<b>Cancer type (CTAB):</b> ", ctab, "</p>",
#         
#         "<p><b>Signature side:</b><br>",
#         "Nomenclature_sig: <b>", nom_s, "</b><br>",
#         "Molecular class: <b>", class_s, "</b><br>",
#         "Signatures: ", sig, "</p>",
#         
#         "<p><b>Interaction side:</b><br>",
#         "Nomenclature_int: <b>", nom_i, "</b><br>",
#         "Molecular class: <b>", class_i, "</b><br>",
#         "Interaction: ", inter, "</p>",
#         
#         "<p><b>Metabolic context:</b><br>",
#         "Metabolism: <b>", metab, "</b><br>",
#         "Pathways: <b>", path, "</b><br>",
#         "Metabolic cell death: <b>", mcd, "</b></p>",
#         
#         "<p><b>Barycenter distance:</b> ",
#         if (!is.na(bary)) sprintf("%.3f", as.numeric(bary)) else "NA",
#         " (", dist_imp, ")</p>",
#         
#         "<p><b>Convex-hull volumes:</b><br>",
#         "sig_hull_vol: ", sig_vol, "<br>",
#         "int_hull_vol: ", int_vol, "<br>",
#         "vol_ratio (sig/int): ",
#         if (!is.na(vr)) sprintf("%.3f", as.numeric(vr)) else "NA",
#         " (", vol_imp, ")</p>"
#       ))
#     })
#     
#     # ---------------- Plot 3D -----------------------------------------
#     output$poly_plot <- renderPlotly({
#       row <- selected_row()
#       req(nrow(row) == 1)
#       
#       # index global do tensor correspondente a esse Circuitries_id
#       m_all <- meta()
#       req("Circuitries_id" %in% names(m_all))
#       idx_global <- match(row$Circuitries_id, m_all$Circuitries_id)
#       req(!is.na(idx_global))
#       
#       # build_circuitry_polytope() e plot_circuitry_polytope()
#       # devem estar definidos no seu código (como no paper).
#       poly_obj <- build_circuitry_polytope(tensor, embedding, index = idx_global)
#       plot_circuitry_polytope(poly_obj)
#     })
#     
#     # ---------------- Download HTML com o plot 3D ----------------------
#     output$dl_poly_html <- downloadHandler(
#       filename = function() {
#         cid <- selected_row()$Circuitries_id
#         paste0("polytope_", gsub("[^A-Za-z0-9]+", "_", cid), ".html")
#       },
#       content = function(file) {
#         row <- selected_row()
#         m_all <- meta()
#         idx_global <- match(row$Circuitries_id, m_all$Circuitries_id)
#         
#         poly_obj <- build_circuitry_polytope(tensor, embedding, index = idx_global)
#         p <- plot_circuitry_polytope(poly_obj)
#         
#         htmlwidgets::saveWidget(as_widget(p), file, selfcontained = TRUE)
#       }
#     )
#   })
# }


# # ============================================================
# # modules/mod_circuitry_polytope.R
# # Módulo completo 3D: filtros → seleção → polytope
# # ============================================================
# 
# mod_circuitry_polytope_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     fluidRow(
#       box(
#         width = 3,
#         title = "Filter regulatory circuitries",
#         status = "info",
#         solidHeader = TRUE,
#         collapsible = TRUE,
#         
#         # FILTROS DO PAPER
#         selectInput(ns("fil_ctab"), "Cancer Type (CTAB):", choices = NULL),
#         selectInput(ns("fil_metabolism"), "Metabolic Superfamily:", choices = NULL),
#         selectInput(ns("fil_pathways"), "Pathways:", choices = NULL),
#         selectInput(ns("fil_mcd"), "Metabolic Cell Death:", choices = NULL),
#         selectInput(ns("fil_dist_imp"), "Distance Regime:", choices = NULL),
#         selectInput(ns("fil_vol_imp"), "Volume Regime:", choices = NULL),
#         
#         hr(),
#         uiOutput(ns("circuitry_selector")),
#         hr(),
#         downloadButton(ns("dl_poly_html"), "Download 3D Polytope (.html)")
#       ),
#       
#       box(
#         width = 9,
#         title = "3D Regulatory Circuitry Polytope",
#         status = "primary",
#         solidHeader = TRUE,
#         plotlyOutput(ns("poly_plot"), height = "600px"),
#         uiOutput(ns("summary_ui"))
#       )
#     )
#   )
# }
# 
# 
# mod_circuitry_polytope_server <- function(id, tensor, embedding) {
#   moduleServer(id, function(input, output, session) {
#     
#     ns <- session$ns
#     
#     meta <- reactive({
#       validate(need(!is.null(tensor), "Tensor not loaded"))
#       validate(need(!is.null(tensor$meta), "Tensor meta missing"))
#       tensor$meta
#     })
#     
#     # -------------------------------------------------------
#     # Inicializa filtros quando meta() é carregado
#     # -------------------------------------------------------
#     observeEvent(meta(), {
#       m <- meta()
#       
#       updateSelectInput(session, "fil_ctab",
#                         choices = sort(unique(m$CTAB)))
#       
#       updateSelectInput(session, "fil_metabolism",
#                         choices = c("Any", sort(unique(m$Metabolism))))
#       
#       updateSelectInput(session, "fil_pathways",
#                         choices = c("Any", sort(unique(m$Pathways))))
#       
#       updateSelectInput(session, "fil_mcd",
#                         choices = c("Any", sort(unique(m$Metabolic_cell_death))))
#       
#       updateSelectInput(session, "fil_dist_imp",
#                         choices = c("Any", sort(unique(m$distance_implication))))
#       
#       updateSelectInput(session, "fil_vol_imp",
#                         choices = c("Any", sort(unique(m$vol_implication))))
#     })
#     
#     # -------------------------------------------------------
#     # Filtragem baseada nos selectInput
#     # -------------------------------------------------------
#     filtered_meta <- reactive({
#       m <- meta()
#       req(nrow(m) > 0)
#       
#       m <- m[m$CTAB == input$fil_ctab, , drop = FALSE]
#       
#       if (input$fil_metabolism != "Any")
#         m <- m[m$Metabolism == input$fil_metabolism,]
#       
#       if (input$fil_pathways != "Any")
#         m <- m[m$Pathways == input$fil_pathways,]
#       
#       if (input$fil_mcd != "Any")
#         m <- m[m$Metabolic_cell_death == input$fil_mcd,]
#       
#       if (input$fil_dist_imp != "Any")
#         m <- m[m$distance_implication == input$fil_dist_imp,]
#       
#       if (input$fil_vol_imp != "Any")
#         m <- m[m$vol_implication == input$fil_vol_imp,]
#       
#       m
#     })
#     
#     # -------------------------------------------------------
#     # UI dinâmico: mostrar somente Circuitries_id válidos
#     # -------------------------------------------------------
#     output$circuitry_selector <- renderUI({
#       m <- filtered_meta()
#       
#       if (nrow(m) == 0)
#         return(p("⚠ No circuitries match filters", style = "color:red;"))
#       
#       selectInput(
#         ns("sel_circ"),
#         "Select Circuitry:",
#         choices = m$Circuitries_id
#       )
#     })
#     
#     selected_index <- reactive({
#       req(filtered_meta(), input$sel_circ)
#       match(input$sel_circ, filtered_meta()$Circuitries_id)
#     })
#     
#     # # -------------------------------------------------------
#     # # RESUMO DO CIRCUITRY ESCOLHIDO
#     # # -------------------------------------------------------
#     # output$summary_ui <- renderUI({
#     #   m <- filtered_meta()
#     #   idx <- selected_index()
#     #   row <- m[idx,]
#     #   
#     #   HTML(paste0(
#     #     "<b>Circuitry:</b> ", row$Circuitries_id, "<br>",
#     #     "<b>Signature:</b> ", row$Nomenclature_sig, "<br>",
#     #     "<b>Interaction:</b> ", row$Nomenclature_int, "<br>",
#     #     "<b>Metabolism:</b> ", row$Metabolism, "<br>",
#     #     "<b>Pathways:</b> ", row$Pathways, "<br>",
#     #     "<b>Cell Death:</b> ", row$Metabolic_cell_death, "<br><br>",
#     #     "<b>Distance regime:</b> ", row$distance_implication, "<br>",
#     #     "<b>Volume regime:</b> ", row$vol_implication
#     #   ))
#     # })
#     
#     # -------------------------------------------------------
#     # RESUMO DO CIRCUITRY ESCOLHIDO (PAR SIGNATURE × INTERACTION)
#     # -------------------------------------------------------
#     output$summary_ui <- renderUI({
#       m   <- filtered_meta()
#       idx <- selected_index()
#       req(m, idx, nrow(m) >= idx)
#       
#       row <- m[idx, , drop = FALSE]
#       
#       # Helper para montar o bloco de uma das metades (signature ou interaction)
#       build_side_block <- function(
#     side_title,
#     nomenclature,
#     molecular_class,
#     components
#       ) {
#         if (is.null(nomenclature) || is.na(nomenclature) || !nzchar(nomenclature)) {
#           return("")
#         }
#         
#         paste0(
#           "<h4>", side_title, "</h4>",
#           "<p>",
#           "<strong>Nomenclature:</strong> ", nomenclature, "<br>",
#           if (!is.null(molecular_class) && !is.na(molecular_class) && nzchar(molecular_class))
#             paste0("<strong>Molecular class:</strong> ", molecular_class, "<br>")
#           else "",
#           if (!is.null(components) && !is.na(components) && nzchar(components))
#             paste0("<strong>Components:</strong> ", components, "<br>")
#           else "",
#           "</p>"
#         )
#       }
#       
#       # Bloco da assinatura (lado da signature)
#       sig_block <- build_side_block(
#         side_title     = "Signature side",
#         nomenclature   = row$Nomenclature_sig,
#         molecular_class = row$Molecular_class_sig,
#         components     = row$Signatures
#       )
#       
#       # Bloco da interação (lado da interaction/regulador)
#       int_block <- build_side_block(
#         side_title     = "Interaction side",
#         nomenclature   = row$Nomenclature_int,
#         molecular_class = row$Molecular_class_int,
#         components     = row$Interaction
#       )
#       
#       # Bloco metabólico/geométrico do par
#       # usa apenas colunas que você já listou que existem
#       cancer_type <- if ("CTAB" %in% names(row)) as.character(row$CTAB) else NA
#       
#       bary_txt <- if ("barycenter_distance" %in% names(row)) {
#         bd <- suppressWarnings(as.numeric(row$barycenter_distance))
#         if (!is.na(bd)) sprintf("%.3f", bd) else "NA"
#       } else {
#         NA
#       }
#       
#       vol_ratio_txt <- if ("vol_ratio" %in% names(row)) {
#         vr <- suppressWarnings(as.numeric(row$vol_ratio))
#         if (!is.na(vr)) sprintf("%.3f", vr) else "NA"
#       } else {
#         NA
#       }
#       
#       conc_summary <- if ("Final_concordance_summary" %in% names(row)) {
#         as.character(row$Final_concordance_summary)
#       } else {
#         NA
#       }
#       
#       meta_block <- paste0(
#         "<h4>Metabolic & geometric context</h4>",
#         "<p>",
#         if (!is.na(cancer_type) && nzchar(cancer_type))
#           paste0("<strong>Cancer type (CTAB):</strong> ", cancer_type, "<br>")
#         else "",
#         
#         if (!is.null(row$Metabolism) && !is.na(row$Metabolism) && nzchar(row$Metabolism))
#           paste0("<strong>Metabolic process:</strong> ", row$Metabolism, "<br>")
#         else "",
#         
#         if (!is.null(row$Pathways) && !is.na(row$Pathways) && nzchar(row$Pathways))
#           paste0("<strong>Metabolic pathway:</strong> ", row$Pathways, "<br>")
#         else "",
#         
#         if (!is.null(row$Metabolic_cell_death) &&
#             !is.na(row$Metabolic_cell_death) &&
#             tolower(row$Metabolic_cell_death) != "unrelated" &&
#             nzchar(row$Metabolic_cell_death))
#           paste0("<strong>Metabolic cell death:</strong> ", row$Metabolic_cell_death, "<br>")
#         else "",
#         
#         if (!is.na(conc_summary) && nzchar(conc_summary))
#           paste0("<strong>Concordance summary:</strong> ", conc_summary, "<br>")
#         else "",
#         
#         if (!is.na(bary_txt))
#           paste0("<strong>Barycenter distance (sig × int):</strong> ", bary_txt, "<br>")
#         else "",
#         
#         if ("distance_implication" %in% names(row) &&
#             !is.na(row$distance_implication) && nzchar(row$distance_implication))
#           paste0("<strong>Distance regime:</strong> ", row$distance_implication, "<br>")
#         else "",
#         
#         if (!is.na(vol_ratio_txt))
#           paste0("<strong>Volume ratio (sig/int):</strong> ", vol_ratio_txt, "<br>")
#         else "",
#         
#         if ("vol_implication" %in% names(row) &&
#             !is.na(row$vol_implication) && nzchar(row$vol_implication))
#           paste0("<strong>Volume regime:</strong> ", row$vol_implication, "<br>")
#         else "",
#         
#         "</p>"
#       )
#       
#       HTML(paste0(
#         "<p><strong>Circuitry ID:</strong> ", row$Circuitries_id, "</p>",
#         sig_block,
#         "<hr style='margin:8px 0;'>",
#         int_block,
#         "<hr style='margin:8px 0;'>",
#         meta_block
#       ))
#     })
#     
#     # -------------------------------------------------------
#     # PLOT 3D
#     # -------------------------------------------------------
#     output$poly_plot <- renderPlotly({
#       m <- filtered_meta()
#       idx <- selected_index()
#       
#       global_idx <- match(m$Circuitries_id[idx], tensor$meta$Circuitries_id)
#       
#       poly <- build_circuitry_polytope(tensor, embedding, index = global_idx)
#       plot_circuitry_polytope(poly)
#     })
#     
#     # -------------------------------------------------------
#     # DOWNLOAD HTML
#     # -------------------------------------------------------
#     output$dl_poly_html <- downloadHandler(
#       filename = function() {
#         paste0("polytope_", input$sel_circ, ".html")
#       },
#       content = function(file) {
#         m <- filtered_meta()
#         idx <- selected_index()
#         global_idx <- match(m$Circuitries_id[idx], tensor$meta$Circuitries_id)
#         
#         poly <- build_circuitry_polytope(tensor, embedding, index = global_idx)
#         p <- plot_circuitry_polytope(poly)
#         
#         htmlwidgets::saveWidget(as_widget(p), file, selfcontained = TRUE)
#       }
#     )
#   })
# }

# # ============================================================
# # modules/mod_circuitry_polytope.R
# # Módulo completo 3D: filtros → seleção → polytope
# # ============================================================
# 
# mod_circuitry_polytope_ui <- function(id) {
#   ns <- NS(id)
#   
#   tagList(
#     fluidRow(
#       box(
#         width = 3,
#         title = "Filter regulatory circuitries",
#         status = "info",
#         solidHeader = TRUE,
#         collapsible = TRUE,
#         
#         # FILTROS DO PAPER
#         selectInput(ns("fil_ctab"), "Cancer Type (CTAB):", choices = NULL),
#         selectInput(ns("fil_metabolism"), "Metabolic Superfamily:", choices = NULL),
#         selectInput(ns("fil_pathways"), "Pathways:", choices = NULL),
#         selectInput(ns("fil_mcd"), "Metabolic Cell Death:", choices = NULL),
#         selectInput(ns("fil_dist_imp"), "Distance Regime:", choices = NULL),
#         selectInput(ns("fil_vol_imp"), "Volume Regime:", choices = NULL),
#         
#         hr(),
#         uiOutput(ns("circuitry_selector")),
#         hr(),
#         downloadButton(ns("dl_poly_html"), "Download 3D Polytope (.html)")
#       ),
#       
#       box(
#         width = 9,
#         title = "3D Regulatory Circuitry Polytope",
#         status = "primary",
#         solidHeader = TRUE,
#         plotlyOutput(ns("poly_plot"), height = "600px"),
#         uiOutput(ns("summary_ui"))
#       )
#     )
#   )
# }
# 
# 
# mod_circuitry_polytope_server <- function(id, tensor, embedding) {
#   moduleServer(id, function(input, output, session) {
#     
#     ns <- session$ns
#     
#     meta <- reactive({
#       validate(need(!is.null(tensor), "Tensor not loaded"))
#       validate(need(!is.null(tensor$meta), "Tensor meta missing"))
#       tensor$meta
#     })
#     
#     # -------------------------------------------------------
#     # Inicializa filtros quando meta() é carregado
#     # -------------------------------------------------------
#     observeEvent(meta(), {
#       m <- meta()
#       
#       updateSelectInput(session, "fil_ctab",
#                         choices = sort(unique(m$CTAB)))
#       
#       updateSelectInput(session, "fil_metabolism",
#                         choices = c("Any", sort(unique(m$Metabolism))))
#       
#       updateSelectInput(session, "fil_pathways",
#                         choices = c("Any", sort(unique(m$Pathways))))
#       
#       updateSelectInput(session, "fil_mcd",
#                         choices = c("Any", sort(unique(m$Metabolic_cell_death))))
#       
#       updateSelectInput(session, "fil_dist_imp",
#                         choices = c("Any", sort(unique(m$distance_implication))))
#       
#       updateSelectInput(session, "fil_vol_imp",
#                         choices = c("Any", sort(unique(m$vol_implication))))
#     })
#     
#     # -------------------------------------------------------
#     # Filtragem baseada nos selectInput
#     # -------------------------------------------------------
#     filtered_meta <- reactive({
#       m <- meta()
#       req(nrow(m) > 0)
#       
#       m <- m[m$CTAB == input$fil_ctab, , drop = FALSE]
#       
#       if (input$fil_metabolism != "Any")
#         m <- m[m$Metabolism == input$fil_metabolism,]
#       
#       if (input$fil_pathways != "Any")
#         m <- m[m$Pathways == input$fil_pathways,]
#       
#       if (input$fil_mcd != "Any")
#         m <- m[m$Metabolic_cell_death == input$fil_mcd,]
#       
#       if (input$fil_dist_imp != "Any")
#         m <- m[m$distance_implication == input$fil_dist_imp,]
#       
#       if (input$fil_vol_imp != "Any")
#         m <- m[m$vol_implication == input$fil_vol_imp,]
#       
#       m
#     })
#     
#     # -------------------------------------------------------
#     # UI dinâmico: mostrar somente Circuitries_id válidos
#     # -------------------------------------------------------
#     output$circuitry_selector <- renderUI({
#       m <- filtered_meta()
#       
#       if (nrow(m) == 0)
#         return(p("⚠ No circuitries match filters", style = "color:red;"))
#       
#       selectInput(
#         ns("sel_circ"),
#         "Select Circuitry:",
#         choices = m$Circuitries_id
#       )
#     })
#     
#     selected_index <- reactive({
#       req(filtered_meta(), input$sel_circ)
#       match(input$sel_circ, filtered_meta()$Circuitries_id)
#     })
#     
#     # # -------------------------------------------------------
#     # # RESUMO DO CIRCUITRY ESCOLHIDO
#     # # -------------------------------------------------------
#     # output$summary_ui <- renderUI({
#     #   m <- filtered_meta()
#     #   idx <- selected_index()
#     #   row <- m[idx,]
#     #   
#     #   HTML(paste0(
#     #     "<b>Circuitry:</b> ", row$Circuitries_id, "<br>",
#     #     "<b>Signature:</b> ", row$Nomenclature_sig, "<br>",
#     #     "<b>Interaction:</b> ", row$Nomenclature_int, "<br>",
#     #     "<b>Metabolism:</b> ", row$Metabolism, "<br>",
#     #     "<b>Pathways:</b> ", row$Pathways, "<br>",
#     #     "<b>Cell Death:</b> ", row$Metabolic_cell_death, "<br><br>",
#     #     "<b>Distance regime:</b> ", row$distance_implication, "<br>",
#     #     "<b>Volume regime:</b> ", row$vol_implication
#     #   ))
#     # })
#     
#     # -------------------------------------------------------
#     # RESUMO DO CIRCUITRY ESCOLHIDO (PAR SIGNATURE × INTERACTION)
#     # -------------------------------------------------------
#     output$summary_ui <- renderUI({
#       m   <- filtered_meta()
#       idx <- selected_index()
#       req(m, idx, nrow(m) >= idx)
#       
#       row <- m[idx, , drop = FALSE]
#       
#       # Helper para montar o bloco de uma das metades (signature ou interaction)
#       build_side_block <- function(
#     side_title,
#     nomenclature,
#     molecular_class,
#     components
#       ) {
#         if (is.null(nomenclature) || is.na(nomenclature) || !nzchar(nomenclature)) {
#           return("")
#         }
#         
#         paste0(
#           "<h4>", side_title, "</h4>",
#           "<p>",
#           "<strong>Nomenclature:</strong> ", nomenclature, "<br>",
#           if (!is.null(molecular_class) && !is.na(molecular_class) && nzchar(molecular_class))
#             paste0("<strong>Molecular class:</strong> ", molecular_class, "<br>")
#           else "",
#           if (!is.null(components) && !is.na(components) && nzchar(components))
#             paste0("<strong>Components:</strong> ", components, "<br>")
#           else "",
#           "</p>"
#         )
#       }
#       
#       # Bloco da assinatura (lado da signature)
#       sig_block <- build_side_block(
#         side_title     = "Signature side",
#         nomenclature   = row$Nomenclature_sig,
#         molecular_class = row$Molecular_class_sig,
#         components     = row$Signatures
#       )
#       
#       # Bloco da interação (lado da interaction/regulador)
#       int_block <- build_side_block(
#         side_title     = "Interaction side",
#         nomenclature   = row$Nomenclature_int,
#         molecular_class = row$Molecular_class_int,
#         components     = row$Interaction
#       )
#       
#       # Bloco metabólico/geométrico do par
#       # usa apenas colunas que você já listou que existem
#       cancer_type <- if ("CTAB" %in% names(row)) as.character(row$CTAB) else NA
#       
#       bary_txt <- if ("barycenter_distance" %in% names(row)) {
#         bd <- suppressWarnings(as.numeric(row$barycenter_distance))
#         if (!is.na(bd)) sprintf("%.3f", bd) else "NA"
#       } else {
#         NA
#       }
#       
#       vol_ratio_txt <- if ("vol_ratio" %in% names(row)) {
#         vr <- suppressWarnings(as.numeric(row$vol_ratio))
#         if (!is.na(vr)) sprintf("%.3f", vr) else "NA"
#       } else {
#         NA
#       }
#       
#       conc_summary <- if ("Final_concordance_summary" %in% names(row)) {
#         as.character(row$Final_concordance_summary)
#       } else {
#         NA
#       }
#       
#       meta_block <- paste0(
#         "<h4>Metabolic & geometric context</h4>",
#         "<p>",
#         if (!is.na(cancer_type) && nzchar(cancer_type))
#           paste0("<strong>Cancer type (CTAB):</strong> ", cancer_type, "<br>")
#         else "",
#         
#         if (!is.null(row$Metabolism) && !is.na(row$Metabolism) && nzchar(row$Metabolism))
#           paste0("<strong>Metabolic process:</strong> ", row$Metabolism, "<br>")
#         else "",
#         
#         if (!is.null(row$Pathways) && !is.na(row$Pathways) && nzchar(row$Pathways))
#           paste0("<strong>Metabolic pathway:</strong> ", row$Pathways, "<br>")
#         else "",
#         
#         if (!is.null(row$Metabolic_cell_death) &&
#             !is.na(row$Metabolic_cell_death) &&
#             tolower(row$Metabolic_cell_death) != "unrelated" &&
#             nzchar(row$Metabolic_cell_death))
#           paste0("<strong>Metabolic cell death:</strong> ", row$Metabolic_cell_death, "<br>")
#         else "",
#         
#         if (!is.na(conc_summary) && nzchar(conc_summary))
#           paste0("<strong>Concordance summary:</strong> ", conc_summary, "<br>")
#         else "",
#         
#         if (!is.na(bary_txt))
#           paste0("<strong>Barycenter distance (sig × int):</strong> ", bary_txt, "<br>")
#         else "",
#         
#         if ("distance_implication" %in% names(row) &&
#             !is.na(row$distance_implication) && nzchar(row$distance_implication))
#           paste0("<strong>Distance regime:</strong> ", row$distance_implication, "<br>")
#         else "",
#         
#         if (!is.na(vol_ratio_txt))
#           paste0("<strong>Volume ratio (sig/int):</strong> ", vol_ratio_txt, "<br>")
#         else "",
#         
#         if ("vol_implication" %in% names(row) &&
#             !is.na(row$vol_implication) && nzchar(row$vol_implication))
#           paste0("<strong>Volume regime:</strong> ", row$vol_implication, "<br>")
#         else "",
#         
#         "</p>"
#       )
#       
#       HTML(paste0(
#         "<p><strong>Circuitry ID:</strong> ", row$Circuitries_id, "</p>",
#         sig_block,
#         "<hr style='margin:8px 0;'>",
#         int_block,
#         "<hr style='margin:8px 0;'>",
#         meta_block
#       ))
#     })
#     
#     # -------------------------------------------------------
#     # PLOT 3D
#     # -------------------------------------------------------
#     output$poly_plot <- renderPlotly({
#       m <- filtered_meta()
#       idx <- selected_index()
#       
#       global_idx <- match(m$Circuitries_id[idx], tensor$meta$Circuitries_id)
#       
#       poly <- build_circuitry_polytope(tensor, embedding, index = global_idx)
#       plot_circuitry_polytope(poly)
#     })
#     
#     # -------------------------------------------------------
#     # DOWNLOAD HTML
#     # -------------------------------------------------------
#     output$dl_poly_html <- downloadHandler(
#       filename = function() {
#         paste0("polytope_", input$sel_circ, ".html")
#       },
#       content = function(file) {
#         m <- filtered_meta()
#         idx <- selected_index()
#         global_idx <- match(m$Circuitries_id[idx], tensor$meta$Circuitries_id)
#         
#         poly <- build_circuitry_polytope(tensor, embedding, index = global_idx)
#         p <- plot_circuitry_polytope(poly)
#         
#         htmlwidgets::saveWidget(as_widget(p), file, selfcontained = TRUE)
#       }
#     )
#   })
# }

# ============================================================
# modules/mod_circuitry_polytope.R
# Módulo completo 3D: filtros → seleção → polytope
# ============================================================

mod_circuitry_polytope_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(
        width = 3,
        title = "Filter regulatory circuitries",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        # FILTROS DO PAPER
        selectInput(ns("fil_ctab"), "Cancer Type (CTAB):", choices = NULL),
        selectInput(ns("fil_metabolism"), "Metabolic Superfamily:", choices = NULL),
        selectInput(ns("fil_pathways"), "Pathways:", choices = NULL),
        selectInput(ns("fil_mcd"), "Metabolic Cell Death:", choices = NULL),
        selectInput(ns("fil_dist_imp"), "Distance Regime:", choices = NULL),
        selectInput(ns("fil_vol_imp"), "Volume Regime:", choices = NULL),
        
        hr(),
        uiOutput(ns("circuitry_selector")),
        hr(),
        downloadButton(ns("dl_poly_html"), "Download 3D Polytope (.html)")
      ),
      
      box(
        width = 9,
        title = "3D Regulatory Circuitry Polytope",
        status = "primary",
        solidHeader = TRUE,
        plotlyOutput(ns("poly_plot"), height = "600px"),
        uiOutput(ns("summary_ui"))
      )
    )
  )
}


mod_circuitry_polytope_server <- function(id, tensor, embedding) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    meta <- reactive({
      validate(need(!is.null(tensor), "Tensor not loaded"))
      validate(need(!is.null(tensor$meta), "Tensor meta missing"))
      tensor$meta
    })
    
    # -------------------------------------------------------
    # Inicializa filtros quando meta() é carregado
    # (once=TRUE reduz recomputação em ambientes limitados)
    # -------------------------------------------------------
    observeEvent(meta(), {
      m <- meta()
      
      updateSelectInput(session, "fil_ctab",
                        choices = sort(unique(m$CTAB)))
      
      updateSelectInput(session, "fil_metabolism",
                        choices = c("Any", sort(unique(m$Metabolism))))
      
      updateSelectInput(session, "fil_pathways",
                        choices = c("Any", sort(unique(m$Pathways))))
      
      updateSelectInput(session, "fil_mcd",
                        choices = c("Any", sort(unique(m$Metabolic_cell_death))))
      
      updateSelectInput(session, "fil_dist_imp",
                        choices = c("Any", sort(unique(m$distance_implication))))
      
      updateSelectInput(session, "fil_vol_imp",
                        choices = c("Any", sort(unique(m$vol_implication))))
    }, once = TRUE)
    
    # -------------------------------------------------------
    # Filtragem baseada nos selectInput
    # -------------------------------------------------------
    filtered_meta <- reactive({
      m <- meta()
      req(nrow(m) > 0, input$fil_ctab)
      
      m <- m[m$CTAB == input$fil_ctab, , drop = FALSE]
      
      if (!is.null(input$fil_metabolism) && input$fil_metabolism != "Any")
        m <- m[m$Metabolism == input$fil_metabolism, , drop = FALSE]
      
      if (!is.null(input$fil_pathways) && input$fil_pathways != "Any")
        m <- m[m$Pathways == input$fil_pathways, , drop = FALSE]
      
      if (!is.null(input$fil_mcd) && input$fil_mcd != "Any")
        m <- m[m$Metabolic_cell_death == input$fil_mcd, , drop = FALSE]
      
      if (!is.null(input$fil_dist_imp) && input$fil_dist_imp != "Any")
        m <- m[m$distance_implication == input$fil_dist_imp, , drop = FALSE]
      
      if (!is.null(input$fil_vol_imp) && input$fil_vol_imp != "Any")
        m <- m[m$vol_implication == input$fil_vol_imp, , drop = FALSE]
      
      m
    })
    
    # IDs filtrados (mais leve do que passar data.frame inteiro para UI)
    filtered_ids <- reactive({
      m <- filtered_meta()
      if (nrow(m) == 0) character(0) else as.character(m$Circuitries_id)
    })
    
    # -------------------------------------------------------
    # UI dinâmico: seletor (usa selectize e carrega choices via updateSelectizeInput server-side)
    # -------------------------------------------------------
    output$circuitry_selector <- renderUI({
      # Mantém o texto exatamente igual
      selectizeInput(
        ns("sel_circ"),
        "Select Circuitry:",
        choices = NULL
      )
    })
    
    # Atualiza choices do selectize de forma server-side (mais leve)
    observeEvent(filtered_ids(), {
      ids <- filtered_ids()
      
      if (length(ids) == 0) {
        updateSelectizeInput(session, "sel_circ",
                             choices = character(0),
                             selected = character(0),
                             server = TRUE)
      } else {
        # mantém seleção se ainda existir; senão, seleciona o primeiro
        current <- isolate(input$sel_circ)
        selected <- if (!is.null(current) && nzchar(current) && current %in% ids) current else ids[1]
        
        updateSelectizeInput(session, "sel_circ",
                             choices = ids,
                             selected = selected,
                             server = TRUE)
      }
    }, ignoreInit = TRUE)
    
    # Mensagem de "sem resultados" sem mudar o texto original
    # (mantém o aviso vermelho quando não há circuitries)
    observe({
      ids <- filtered_ids()
      if (length(ids) == 0) {
        # substitui o UI do selectize por um aviso (igual ao seu texto)
        output$circuitry_selector <- renderUI({
          p("⚠ No circuitries match filters", style = "color:red;")
        })
      } else {
        # restaura o selectizeInput quando voltar a ter opções
        output$circuitry_selector <- renderUI({
          selectizeInput(
            ns("sel_circ"),
            "Select Circuitry:",
            choices = NULL
          )
        })
      }
    })
    
    # -------------------------------------------------------
    # Índice global direto do tensor (evita match duplicado)
    # -------------------------------------------------------
    global_idx <- reactive({
      req(input$sel_circ)
      match(input$sel_circ, tensor$meta$Circuitries_id)
    })
    
    # -------------------------------------------------------
    # RESUMO DO CIRCUITRY ESCOLHIDO (PAR SIGNATURE × INTERACTION)
    # (o seu código original, sem alterações de texto)
    # -------------------------------------------------------
    output$summary_ui <- renderUI({
      m <- filtered_meta()
      req(nrow(m) > 0, input$sel_circ)
      
      idx <- match(input$sel_circ, m$Circuitries_id)
      req(!is.na(idx))
      
      row <- m[idx, , drop = FALSE]
      
      build_side_block <- function(
    side_title,
    nomenclature,
    molecular_class,
    components
      ) {
        if (is.null(nomenclature) || is.na(nomenclature) || !nzchar(nomenclature)) {
          return("")
        }
        
        paste0(
          "<h4>", side_title, "</h4>",
          "<p>",
          "<strong>Nomenclature:</strong> ", nomenclature, "<br>",
          if (!is.null(molecular_class) && !is.na(molecular_class) && nzchar(molecular_class))
            paste0("<strong>Molecular class:</strong> ", molecular_class, "<br>")
          else "",
          if (!is.null(components) && !is.na(components) && nzchar(components))
            paste0("<strong>Components:</strong> ", components, "<br>")
          else "",
          "</p>"
        )
      }
      
      sig_block <- build_side_block(
        side_title      = "Signature side",
        nomenclature    = row$Nomenclature_sig,
        molecular_class = row$Molecular_class_sig,
        components      = row$Signatures
      )
      
      int_block <- build_side_block(
        side_title      = "Interaction side",
        nomenclature    = row$Nomenclature_int,
        molecular_class = row$Molecular_class_int,
        components      = row$Interaction
      )
      
      cancer_type <- if ("CTAB" %in% names(row)) as.character(row$CTAB) else NA
      
      bary_txt <- if ("barycenter_distance" %in% names(row)) {
        bd <- suppressWarnings(as.numeric(row$barycenter_distance))
        if (!is.na(bd)) sprintf("%.3f", bd) else "NA"
      } else {
        NA
      }
      
      vol_ratio_txt <- if ("vol_ratio" %in% names(row)) {
        vr <- suppressWarnings(as.numeric(row$vol_ratio))
        if (!is.na(vr)) sprintf("%.3f", vr) else "NA"
      } else {
        NA
      }
      
      conc_summary <- if ("Final_concordance_summary" %in% names(row)) {
        as.character(row$Final_concordance_summary)
      } else {
        NA
      }
      
      meta_block <- paste0(
        "<h4>Metabolic & geometric context</h4>",
        "<p>",
        if (!is.na(cancer_type) && nzchar(cancer_type))
          paste0("<strong>Cancer type (CTAB):</strong> ", cancer_type, "<br>")
        else "",
        
        if (!is.null(row$Metabolism) && !is.na(row$Metabolism) && nzchar(row$Metabolism))
          paste0("<strong>Metabolic process:</strong> ", row$Metabolism, "<br>")
        else "",
        
        if (!is.null(row$Pathways) && !is.na(row$Pathways) && nzchar(row$Pathways))
          paste0("<strong>Metabolic pathway:</strong> ", row$Pathways, "<br>")
        else "",
        
        if (!is.null(row$Metabolic_cell_death) &&
            !is.na(row$Metabolic_cell_death) &&
            tolower(row$Metabolic_cell_death) != "unrelated" &&
            nzchar(row$Metabolic_cell_death))
          paste0("<strong>Metabolic cell death:</strong> ", row$Metabolic_cell_death, "<br>")
        else "",
        
        if (!is.na(conc_summary) && nzchar(conc_summary))
          paste0("<strong>Concordance summary:</strong> ", conc_summary, "<br>")
        else "",
        
        if (!is.na(bary_txt))
          paste0("<strong>Barycenter distance (sig × int):</strong> ", bary_txt, "<br>")
        else "",
        
        if ("distance_implication" %in% names(row) &&
            !is.na(row$distance_implication) && nzchar(row$distance_implication))
          paste0("<strong>Distance regime:</strong> ", row$distance_implication, "<br>")
        else "",
        
        if (!is.na(vol_ratio_txt))
          paste0("<strong>Volume ratio (sig/int):</strong> ", vol_ratio_txt, "<br>")
        else "",
        
        if ("vol_implication" %in% names(row) &&
            !is.na(row$vol_implication) && nzchar(row$vol_implication))
          paste0("<strong>Volume regime:</strong> ", row$vol_implication, "<br>")
        else "",
        
        "</p>"
      )
      
      HTML(paste0(
        "<p><strong>Circuitry ID:</strong> ", row$Circuitries_id, "</p>",
        sig_block,
        "<hr style='margin:8px 0;'>",
        int_block,
        "<hr style='margin:8px 0;'>",
        meta_block
      ))
    })
    
    # -------------------------------------------------------
    # Polytope: calcula só quando sel_circ muda (reduz CPU)
    # -------------------------------------------------------
    poly_obj <- eventReactive(input$sel_circ, {
      idx <- global_idx()
      req(!is.na(idx))
      build_circuitry_polytope(tensor, embedding, index = idx)
    }, ignoreInit = TRUE)
    
    output$poly_plot <- renderPlotly({
      req(poly_obj())
      plot_circuitry_polytope(poly_obj())
    })
    
    # Se o módulo estiver em tab/accordion oculto, evita gastar CPU enquanto escondido
    outputOptions(output, "poly_plot", suspendWhenHidden = TRUE)
    
    # -------------------------------------------------------
    # DOWNLOAD HTML (reusa poly_obj quando disponível)
    # -------------------------------------------------------
    output$dl_poly_html <- downloadHandler(
      filename = function() {
        paste0("polytope_", input$sel_circ, ".html")
      },
      content = function(file) {
        req(input$sel_circ)
        
        poly <- poly_obj()
        if (is.null(poly)) {
          idx <- global_idx()
          req(!is.na(idx))
          poly <- build_circuitry_polytope(tensor, embedding, index = idx)
        }
        
        p <- plot_circuitry_polytope(poly)
        htmlwidgets::saveWidget(as_widget(p), file, selfcontained = TRUE)
      }
    )
    
  })
}

