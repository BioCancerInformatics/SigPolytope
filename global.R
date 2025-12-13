# # =========================================================
# # Load Required Packages
# # =========================================================
# suppressPackageStartupMessages({
#   library(shiny)
#   library(shinycssloaders)
#   library(shinydashboard)
#   
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(purrr)
#   library(readr)
#   library(readxl)
#   library(glue)
#   library(memoise)
#   library(rio)
#   
#   library(ggplot2)
#   library(ggsci)
#   library(cowplot)
#   library(grid)
#   library(gridtext)
#   library(geometry)
#   library(plotly)
#   
#   library(igraph)
#   library(tidygraph)
#   library(ggraph)
#   library(visNetwork)
#   
#   library(survival)
#   library(survminer)
#   library(fmsb)
#   library(ggradar)
#   
#   library(DT)
#   library(htmltools)
#   library(htmlwidgets)
#   library(UCSCXenaShiny)
#   library(UCSCXenaTools)
# })
# 
# options(shiny.maxRequestSize = 500 * 1024^2)
# options(shiny.sanitize.errors = TRUE)
# 
# # ---------------------- Carregamento lazy do All_data ----------------------
# load_all_data <- memoise(function() {
#   file_path <- file.path("data", "All_data.rds")
#   if (!file.exists(file_path)) {
#     stop("Arquivo 'data/All_data.rds' não encontrado.")
#   }
#   message("Carregando All_data.rds...")
#   readRDS(file_path)
# })
# 
# # TCGA types (necessário para módulo de interações)
# tcga_types <- readxl::read_excel("data/TCGA_Cancer_types.xlsx")
# 
# # ---------------------- Importação de Módulos ----------------------
# source("modules/module_search_your_target.R")
# source("modules/module_nomenclature.R")
# source("modules/module_signature.R")
# source("modules/module_all_signatures.R")
# source("modules/module_meaningful_interaction.R")
# source("modules/module_developers.R")
# 
# # Funções de geometria (tensor, PCA, polytope, add_geometry_metadata)
# source("R/geometry_polytope_funs.R")
# 
# # Módulo 3D
# source("modules/module_circuitry_polytope.R")
# 
# # ---------------------- Carregar dataset de circuitries ----------------------
# circuitries <- readRDS("data/Regulatory_circuitries_geometric.rds")
# # se o correto for 'circuitries_final.rds', troque a linha acima
# 
# # ---------------------- Construir tensor 18D ----------------------
# tensor <- build_geometric_tensor_from_circuitries(circuitries)
# 
# # ---------------------- PCA (18D -> 3D) --------------------------
# embedding <- compute_pca_embedding_from_tensor(tensor)



# =========================================================
# global.R (otimizado para CPU/RAM + fallback seguro)
# =========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(shinydashboard)
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(readr)
  library(readxl)
  library(glue)
  library(memoise)
  
  library(htmltools)
  library(htmlwidgets)
  
  library(ggplot2)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(visNetwork)
  library(DT)
  
  library(plotly)
  library(geometry)
})

options(shiny.maxRequestSize = 100 * 1024^2)
options(shiny.sanitize.errors = TRUE)

# ---------------------- All_data.rds ----------------------
load_all_data <- memoise(function() {
  file_path <- file.path("data", "All_data.rds")
  if (!file.exists(file_path)) stop("Arquivo 'data/All_data.rds' não encontrado.")
  message("Carregando All_data.rds...")
  readRDS(file_path)
})

# ---------------------- TCGA types (RDS preferido, XLSX fallback) ----------------------
load_tcga_types <- memoise(function() {
  rds_path  <- file.path("data", "TCGA_Cancer_types.rds")
  xlsx_path <- file.path("data", "TCGA_Cancer_types.xlsx")
  
  if (file.exists(rds_path)) {
    message("Carregando TCGA_Cancer_types.rds...")
    return(readRDS(rds_path))
  }
  
  if (file.exists(xlsx_path)) {
    message("Carregando TCGA_Cancer_types.xlsx (fallback)...")
    return(readxl::read_excel(xlsx_path))
  }
  
  stop("Nenhum arquivo de TCGA types encontrado (nem .rds nem .xlsx em /data).")
})

# ---------------------- Tensor + Embedding (RDS preferido, build fallback) ----------------------
get_tensor_embedding <- memoise(function() {
  te_path <- file.path("data", "tensor_embedding.rds")
  if (file.exists(te_path)) {
    message("Carregando tensor/embedding pré-computados (tensor_embedding.rds)...")
    te <- readRDS(te_path)
    if (!is.list(te) || is.null(te$tensor) || is.null(te$embedding)) {
      stop("'tensor_embedding.rds' deve ser um list(tensor=..., embedding=...).")
    }
    return(te)
  }
  
  # Fallback: construir no runtime (pesado; bom só localmente)
  circuitries_path <- file.path("data", "Regulatory_circuitries_geometric.rds")
  if (!file.exists(circuitries_path)) {
    stop("Arquivo 'data/Regulatory_circuitries_geometric.rds' não encontrado.")
  }
  
  message("Construindo tensor/embedding no runtime (fallback pesado)...")
  circuitries <- readRDS(circuitries_path)
  
  # Essas funções vêm dos seus sources
  tensor <- build_geometric_tensor_from_circuitries(circuitries)
  embedding <- compute_pca_embedding_from_tensor(tensor)
  
  list(tensor = tensor, embedding = embedding)
})

# ---------------------- Utilitário: pré-computar RDS localmente ----------------------
precompute_assets <- function(save_tcga = TRUE, save_te = TRUE) {
  dir.create("data", showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(save_tcga)) {
    xlsx_path <- file.path("data", "TCGA_Cancer_types.xlsx")
    rds_path  <- file.path("data", "TCGA_Cancer_types.rds")
    if (file.exists(xlsx_path) && !file.exists(rds_path)) {
      message("Gerando TCGA_Cancer_types.rds...")
      tc <- readxl::read_excel(xlsx_path)
      saveRDS(tc, rds_path)
    }
  }
  
  if (isTRUE(save_te)) {
    te_path <- file.path("data", "tensor_embedding.rds")
    if (!file.exists(te_path)) {
      message("Gerando tensor_embedding.rds (pode demorar)...")
      te <- get_tensor_embedding()  # isso vai cair no fallback pesado se não existir
      saveRDS(te, te_path)
    }
  }
  
  invisible(TRUE)
}

# ---------------------- Importação de módulos ----------------------
source("modules/module_meaningful_interaction.R")
source("modules/module_developers.R")

source("R/geometry_polytope_funs.R")
source("modules/module_circuitry_polytope.R")

