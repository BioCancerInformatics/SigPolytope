# # =========================================================
# # server.R
# # =========================================================
# 
# # (Opcional) Contador de visitas: arquivo local
# counter_file <- "counter.txt"
# 
# update_counter <- function() {
#   # Em shinyapps.io, o sistema de arquivos é efêmero,
#   # mas é gravável durante a sessão. Ainda assim,
#   # envolvemos em tryCatch para evitar que erros de I/O
#   # impeçam o app de subir.
#   tryCatch({
#     if (!file.exists(counter_file)) {
#       writeLines("0", counter_file)
#     }
#     count <- as.integer(readLines(counter_file))
#     if (is.na(count)) count <- 0
#     count <- count + 1
#     writeLines(as.character(count), counter_file)
#     count
#   }, error = function(e) {
#     message("Falha ao atualizar o contador: ", conditionMessage(e))
#     NA_integer_
#   })
# }
# 
# server <- function(input, output, session) {
#   session$allowReconnect(TRUE)  # ok, pode manter
#   
#   # ---------------------- Carregar all_data ----------------------
#   # Aqui usamos a função memoizada definida no global.R
#   all_data <- NULL
#   all_data <- tryCatch(
#     {
#       load_all_data()
#     },
#     error = function(e) {
#       # Essa mensagem aparecerá no log do shinyapps.io
#       message("Erro ao carregar all_data.rds: ", conditionMessage(e))
#       NULL
#     }
#   )
#   
#   if (is.null(all_data)) {
#     # Se chegou aqui, provavelmente o arquivo está faltando
#     # ou estourando limite. Você pode exibir uma mensagem amigável
#     # na interface se quiser.
#     showModal(
#       modalDialog(
#         title = "Erro de inicialização",
#         "Não foi possível carregar os dados necessários (all_data.rds). 
#         Verifique se o arquivo foi incluído corretamente no deploy 
#         e se o tamanho não excede os limites do shinyapps.io.",
#         easyClose = TRUE
#       )
#     )
#     # Podemos retornar cedo, evitando erros em cascata:
#     return(invisible(NULL))
#   }
# 
#   # # ---------------------- Contador de visitas ----------------------
#   # visitor_count <- update_counter()
#   # 
#   # output$visitor_count <- renderText({
#   #   if (is.na(visitor_count)) {
#   #     "Total: (contador indisponível no servidor)"
#   #   } else {
#   #     paste("Total:", visitor_count)
#   #   }
#   # })
#   
#   # ---------------------- Reativos principais ----------------------
#   # Aqui assumo que all_data é uma lista com esses nomes.
#   # Se o nome de algum elemento for diferente, é só ajustar.
#   
#   search_your_target      <- reactiveVal(all_data$Search_your_target)
#   Target                  <- reactiveVal(all_data$Dataset_S3)
#   meaningful_interaction  <- reactiveVal(all_data$Dataset_S4)
# 
#   
#   # ---------------------- Chamada dos módulos ----------------------
#   # DASHBOARD: Search Your Gene + Regulatory circuitry integrado
#   mod_search_your_target_server(
#     id                 = "home",
#     search_dataset     = search_your_target,
#     Target_dataset     = Target,
#     tcga_types         = tcga_types,
#     regulatory_dataset = meaningful_interaction
#   )
#   
#   
#   
#   # 
#   # output$overview_poly_demo <- plotly::renderPlotly({
#   #   
#   #   # 1) garantir que o tensor e embedding já existem
#   #   req(exists("tensor", inherits = TRUE))
#   #   req(exists("embedding", inherits = TRUE))
#   #   
#   #   tensor_local <- get("tensor", inherits = TRUE)
#   #   embedding_local <- get("embedding", inherits = TRUE)
#   #   
#   #   req(!is.null(tensor_local$meta))
#   #   req(nrow(tensor_local$meta) > 0)
#   #   
#   #   # 2) escolher UM circuitry leve (exemplo automático)
#   #   # use o primeiro para evitar sample pesado no load
#   #   idx <- 1
#   #   
#   #   # 3) construir o polytope
#   #   poly <- build_circuitry_polytope(
#   #     tensor   = tensor_local,
#   #     embedding = embedding_local,
#   #     index    = idx
#   #   )
#   #   
#   #   # 4) plotar
#   #   plot_circuitry_polytope(
#   #     poly,
#   #     show_vertices = TRUE,
#   #     show_hull = TRUE,
#   #     show_barycenters = TRUE,
#   #     title = "Example geometric regulatory circuitry"
#   #   )
#   # })
#   # 
#   # 
#   
#   # tensor e embedding já estão no ambiente global (vieram do global.R)
#   mod_circuitry_polytope_server(
#     id        = "poly3d",
#     tensor    = tensor,       # objeto que tem tensor_obj$meta = seu data.frame meta
#     embedding = embedding      # matriz / data.frame de coordenadas do espaço latente
#   )
#   
#   
#   mod_meaningful_interaction_server("Mean_int", meaningful_interaction)
#   
#   mod_developers_server("Developers")
#   
#   # ---------------------- Limpeza em fim de sessão ----------------------
#   # NÃO usar rm(list = ls()) aqui: isso mexe no ambiente do worker inteiro.
#   # Use apenas o que for realmente necessário (ex.: desconectar banco, fechar conexões, etc.).
#   session$onSessionEnded(function() {
#     message("Session ended.")
#     # Se você criar objetos específicos da sessão (ex. conexões DB),
#     # limpe-os aqui explicitamente.
#     gc()
#   })
# }
# ]


# =========================================================
# server.R (otimizado para CPU/RAM + fallback seguro)
# =========================================================

server <- function(input, output, session) {
  session$allowReconnect(TRUE)
  
  # ---------------------- Carregar all_data ----------------------
  all_data <- tryCatch(
    load_all_data(),
    error = function(e) {
      message("Erro ao carregar All_data.rds: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(all_data)) {
    showModal(
      modalDialog(
        title = "Erro de inicialização",
        "Não foi possível carregar os dados necessários (all_data.rds). 
        Verifique se o arquivo foi incluído corretamente no deploy 
        e se o tamanho não excede os limites do shinyapps.io.",
        easyClose = TRUE
      )
    )
    return(invisible(NULL))
  }
  
  # ---------------------- tcga_types (RDS preferido, XLSX fallback) ----------------------
  tcga_types <- tryCatch(
    load_tcga_types(),
    error = function(e) {
      message("Erro ao carregar TCGA types: ", conditionMessage(e))
      NULL
    }
  )
  
  # ---------------------- Extrair só o que precisa e soltar a lista grande ----------------------
  dataset_S4 <- all_data$Dataset_S4
  rm(all_data)
  
  meaningful_interaction <- reactive(dataset_S4)
  
  # ---------------------- Tensor/Embedding ----------------------
  te <- tryCatch(
    get_tensor_embedding(),
    error = function(e) {
      message("Erro ao obter tensor/embedding: ", conditionMessage(e))
      NULL
    }
  )
  
  # Sobe o resto do app mesmo se o 3D falhar
  if (!is.null(te)) {
    mod_circuitry_polytope_server(
      id        = "poly3d",
      tensor    = te$tensor,
      embedding = te$embedding
    )
  } else {
    showModal(
      modalDialog(
        title = "Erro de inicialização",
        "Não foi possível carregar o tensor/embedding do módulo 3D.",
        easyClose = TRUE
      )
    )
  }
  
  mod_meaningful_interaction_server(
    id         = "Mean_int",
    dataset    = meaningful_interaction,
    tcga_types = tcga_types
  )
  
  mod_developers_server("Developers")
  
  session$onSessionEnded(function() {
    message("Session ended.")
  })
}
