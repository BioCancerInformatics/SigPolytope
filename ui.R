#' ui <- dashboardPage(
#'   skin = "black",
#'   
#'   dashboardHeader(
#'     title = tagList(
#'       span("SigPolytope Shiny", style = "font-weight:600;")
#'     ),
#'     titleWidth = 260
#'   ),
#'   
#'   dashboardSidebar(
#'     width = 260,
#'     sidebarMenu(
#'       id = "main_menu",
#'       menuItem("Overview & Concept", tabName = "overview", icon = icon("home")),
#'       menuItem("Regulatory Circuitries (3D)", tabName = "geometry", icon = icon("cube")),
#'       menuItem("Regulatory Circuitries (network)", tabName = "reg_circ", icon = icon("share-alt")),
#'       menuItem("About & Citation", tabName = "about", icon = icon("info-circle")),
#'       menuItem("GitHub", tabName = "github_link", icon = icon("github"))
#'     )
#'   ),
#'   
#'   dashboardBody(
#'     tags$head(
#'       tags$style(HTML("
#'         .content-wrapper, .right-side { background-color: #f7f7f7; }
#'         .box { border-radius: 10px; }
#'         .small-box { border-radius: 10px !important; }
#' 
#'         /* HERO */
#'         .overview-hero { min-height: 640px; }
#'         .overview-hero .box-body { padding: 26px 22px; }
#' 
#'         .overview-title {
#'           font-size: 28px;
#'           font-weight: 700;
#'           color: #111;
#'           margin-bottom: 10px;
#'         }
#' 
#'         .overview-subtitle {
#'           font-size: 16px;
#'           color: #222;
#'           margin-bottom: 14px;
#'           line-height: 1.45;
#'         }
#' 
#'         .overview-text {
#'           font-size: 14px;
#'           color: #333;
#'           line-height: 1.75;
#'         }
#' 
#'         .plot-note {
#'           margin-top: 8px;
#'           font-size: 12px;
#'           color: #666;
#'         }
#' 
#'         /* Loading overlay */
#'         #app-loading-overlay {
#'           position: fixed; top: 0; left: 0;
#'           width: 100vw; height: 100vh;
#'           background: rgba(255,255,255,0.97);
#'           z-index: 99999;
#'           display: flex;
#'           align-items: center;
#'           justify-content: center;
#'           flex-direction: column;
#'           font-family: Arial, sans-serif;
#'         }
#' 
#'         #app-loading-title {
#'           font-size: 26px;
#'           font-weight: 700;
#'           color: #111;
#'           margin-bottom: 10px;
#'         }
#' 
#'         #app-loading-subtitle {
#'           font-size: 14px;
#'           color: #555;
#'           max-width: 720px;
#'           text-align: center;
#'           line-height: 1.5;
#'           margin-bottom: 18px;
#'         }
#' 
#'         .spinner {
#'           width: 44px; height: 44px;
#'           border: 5px solid #ddd;
#'           border-top: 5px solid #111;
#'           border-radius: 50%;
#'           animation: spin 0.9s linear infinite;
#'         }
#' 
#'         @keyframes spin {
#'           0% { transform: rotate(0deg); }
#'           100% { transform: rotate(360deg); }
#'         }
#'       "))
#'     ),
#'     
#'     # Loading overlay
#'     tags$div(
#'       id = "app-loading-overlay",
#'       tags$div(id = "app-loading-title", "SigPolytope Shiny"),
#'       tags$div(
#'         id = "app-loading-subtitle",
#'         HTML("A geometric atlas of multi-omic regulatory circuitries.")
#'       ),
#'       tags$div(class = "spinner"),
#'       tags$div(style="margin-top:12px;color:#777;font-size:12px;", "Loading…")
#'     ),
#'     
#'     tabItems(
#'       
#'       # ---------------- OVERVIEW ----------------
#'       tabItem(
#'         tabName = "overview",
#'         
#'         # HERO
#'         fluidRow(
#'           box(
#'             width = 12,
#'             status = "primary",
#'             solidHeader = TRUE,
#'             class = "overview-hero",
#'             
#'             fluidRow(
#'               column(
#'                 width = 5,
#'                 tags$div(class = "overview-title", "SigPolytope Shiny"),
#'                 tags$div(
#'                   class = "overview-subtitle",
#'                   "A geometric atlas of OncoMetabolism multi-omic regulatory circuitries."
#'                 ),
#'                 tags$div(
#'                   class = "overview-text",
#'                   "Explore omic signatures as structured geometric objects rather than vectors or gene lists. ",
#'                   "SigPolytope enables interactive inspection of convex-polytope representations that preserve ",
#'                   "multidimensional organization, revealing concordance, discordance, and latent complexity across ",
#'                   "metabolic regulatory circuitries."
#'                 )
#'               ),
#'               column(
#'                 width = 7,
#'                 plotly::plotlyOutput("overview_poly_demo", height = "560px"),
#'                 tags$div(
#'                   class = "plot-note",
#'                   "Automatic example rendered at startup (illustrative preview)."
#'                 )
#'               )
#'             )
#'           )
#'         ),
#'         
#'         # -------- RESTORED BOXES --------
#'         fluidRow(
#'           box(
#'             width = 6,
#'             title = "Quick start",
#'             status = "info",
#'             solidHeader = TRUE,
#'             tags$ol(
#'               style = "font-size:13px; line-height:1.65;",
#'               tags$li("Open Regulatory Circuitries (3D)."),
#'               tags$li("Filter by cancer type (CTAB), metabolism and pathway."),
#'               tags$li("Refine using discordance (barycenter distance) and volume regimes."),
#'               tags$li("Select a Circuitries_id to render the paired polytope."),
#'               tags$li("Use the network tab for interaction topology and concordance.")
#'             )
#'           ),
#'           
#'           box(
#'             width = 6,
#'             title = "What you are seeing",
#'             status = "info",
#'             solidHeader = TRUE,
#'             tags$ul(
#'               style = "font-size:13px; line-height:1.65;",
#'               tags$li("Two geometric objects: signature (sig) and interaction (int)."),
#'               tags$li("3D projection of an 18D latent clinical space."),
#'               tags$li("Convex hulls summarizing latent variability."),
#'               tags$li("Barycenter distance capturing multidimensional discordance."),
#'               tags$li("Hull volumes reflecting complexity and asymmetry.")
#'             )
#'           )
#'         )
#'       ),
#'       
#'       # ---------------- GEOMETRY ----------------
#'       tabItem(
#'         tabName = "geometry",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             status = "primary",
#'             solidHeader = TRUE,
#'             title = "3D Geometric Explorer of Regulatory Circuitries",
#'             mod_circuitry_polytope_ui("poly3d")
#'           )
#'         )
#'       ),
#'       
#'       # ---------------- NETWORK ----------------
#'       tabItem(
#'         tabName = "reg_circ",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "Regulatory circuitries — network view",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             mod_meaningful_interaction_ui("Mean_int")
#'           )
#'         )
#'       ),
#'       
#'       # ---------------- ABOUT ----------------
#'       tabItem(
#'         tabName = "about",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "How to cite this geometric framework",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             p("If you use this geometric representation or any output generated by this dashboard, please cite:"),
#'             p(
#'               "Nogueira, H. A. C.; Souza, E. R.; Lopes, V. S.; Medina-Acosta, E. (2025). ",
#'               em("A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer."),
#'               " Preprint. DOI: ",
#'               a(
#'                 "10.1101/2025.11.15.688631",
#'                 href = "https://doi.org/10.1101/2025.11.15.688631",
#'                 target = "_blank"
#'               )
#'             )
#'           )
#'         )
#'       )
#'     ),
#'     
#'     tags$script(HTML("
#'       $(document).on('shiny:connected', function(){
#'         setTimeout(function(){
#'           $('#app-loading-overlay').fadeOut(600);
#'         }, 10000);
#'       });
#'     "))
#'   )
#' )
#' 


#' # =========================================================
#' # UI - Geometric Multidimensional Representation Dashboard
#' # =========================================================
#' 
#' ui <- dashboardPage(
#'   skin = "black",
#'   
#'   dashboardHeader(
#'     title = tagList(
#'       span("SigPolytope Shiny", style = "font-weight:600;"),
#'       HTML('<span style="font-size:11px; margin-left:4px;">(OncoMetabolismGPS – Geometry)</span>')
#'     ),
#'     titleWidth = 260
#'   ),
#'   
#'   dashboardSidebar(
#'     width = 260,
#'     sidebarMenu(
#'       id = "main_menu",
#'       
#'       menuItem("Overview & Concept", tabName = "overview", icon = icon("home")),
#'       menuItem("Regulatory Circuitries (3D)", tabName = "geometry", icon = icon("cube")),
#'       menuItem("Regulatory Circuitries (network)", tabName = "reg_circ", icon = icon("share-alt")),
#'       menuItem("About & Citation", tabName = "about", icon = icon("info-circle")),
#'       menuItem("GitHub", tabName = "github_link", icon = icon("github"))
#'     )
#'   ),
#'   
#'   dashboardBody(
#'     
#'     # ------------------------------------------------------------------
#'     # (ALTERAÇÃO 1) CSS + JS para: Loading 5s + Splash "SigPolytope"
#'     # ------------------------------------------------------------------
#'     tags$head(
#'       tags$style(HTML("
#'         .content-wrapper, .right-side { background-color: #f7f7f7; }
#'         .box { border-radius: 10px; }
#'         .small-box { border-radius: 10px !important; }
#' 
#'         /* =========================
#'            Loading overlay (5s)
#'            ========================= */
#'         #app-loading-overlay {
#'           position: fixed; top: 0; left: 0;
#'           width: 100vw; height: 100vh;
#'           background: rgba(255,255,255,0.97);
#'           z-index: 99999;
#'           display: flex;
#'           align-items: center;
#'           justify-content: center;
#'           flex-direction: column;
#'           font-family: Arial, sans-serif;
#'         }
#' 
#'         #app-loading-title {
#'           font-size: 26px;
#'           font-weight: 700;
#'           color: #111;
#'           margin-bottom: 10px;
#'         }
#' 
#'         #app-loading-subtitle {
#'           font-size: 14px;
#'           color: #555;
#'           max-width: 720px;
#'           text-align: center;
#'           line-height: 1.5;
#'           margin-bottom: 18px;
#'         }
#' 
#'         .spinner {
#'           width: 44px; height: 44px;
#'           border: 5px solid #ddd;
#'           border-top: 5px solid #111;
#'           border-radius: 50%;
#'           animation: spin 0.9s linear infinite;
#'         }
#' 
#'         @keyframes spin {
#'           0% { transform: rotate(0deg); }
#'           100% { transform: rotate(360deg); }
#'         }
#' 
#'         /* =========================
#'            Splash overlay (click-to-dismiss)
#'            ========================= */
#'         #sigpolytope_splash {
#'           position: fixed;
#'           inset: 0;
#'           z-index: 99998; /* abaixo do loading (99999) */
#'           background: rgba(0,0,0,0.70);
#'           display: none;  /* aparece após loading */
#'           align-items: center;
#'           justify-content: center;
#'           padding: 24px;
#'         }
#' 
#'         #sigpolytope_splash .card {
#'           max-width: 980px;
#'           width: 100%;
#'           background: #ffffff;
#'           border-radius: 14px;
#'           padding: 22px 24px;
#'           box-shadow: 0 10px 30px rgba(0,0,0,0.25);
#'         }
#' 
#'         #sigpolytope_splash h2 {
#'           margin-top: 0;
#'           font-weight: 700;
#'         }
#' 
#'         #sigpolytope_splash p {
#'           font-size: 15px;
#'           line-height: 1.45;
#'           margin-bottom: 10px;
#'           color: #222;
#'         }
#' 
#'         #sigpolytope_splash .hint {
#'           margin-top: 14px;
#'           font-size: 13px;
#'           color: #666;
#'         }
#'       ")),
#'       tags$script(HTML("
#'         (function(){
#' 
#'           function hideSplash(){
#'             var el = document.getElementById('sigpolytope_splash');
#'             if(el){ el.style.display = 'none'; }
#'             document.removeEventListener('click', hideSplash, true);
#'             document.removeEventListener('keydown', hideSplash, true);
#'             document.removeEventListener('wheel', hideSplash, true);
#'             document.removeEventListener('touchstart', hideSplash, true);
#'           }
#' 
#'           function showSplash(){
#'             var el = document.getElementById('sigpolytope_splash');
#'             if(el){ el.style.display = 'flex'; }
#'             // some no primeiro clique/interação
#'             document.addEventListener('click', hideSplash, true);
#'             document.addEventListener('keydown', hideSplash, true);
#'             document.addEventListener('wheel', hideSplash, true);
#'             document.addEventListener('touchstart', hideSplash, true);
#'           }
#' 
#'           function hideLoadingAndShowSplash(){
#'             var load = document.getElementById('app-loading-overlay');
#'             if(load){ load.style.display = 'none'; }
#'             showSplash();
#'           }
#' 
#'           // mantém loading por 5s e depois mostra splash
#'           document.addEventListener('DOMContentLoaded', function(){
#'             setTimeout(hideLoadingAndShowSplash, 5000);
#'           });
#' 
#'         })();
#'       "))
#'     ),
#'     
#'     # ------------------------------------------------------------------
#'     # (ALTERAÇÃO 2) Loading overlay (aparece imediatamente)
#'     # ------------------------------------------------------------------
#'     tags$div(
#'       id = "app-loading-overlay",
#'       tags$div(id = "app-loading-title", "SigPolytope Shiny"),
#'       tags$div(
#'         id = "app-loading-subtitle",
#'         HTML("A geometric atlas of multi-omic regulatory circuitries.")
#'       ),
#'       tags$div(class = "spinner"),
#'       tags$div(style="margin-top:12px;color:#777;font-size:12px;", "Loading…")
#'     ),
#'     
#'     # ------------------------------------------------------------------
#'     # (ALTERAÇÃO 3) Splash / mensagem inicial (aparece após 5s; some no clique)
#'     # ------------------------------------------------------------------
#'     tags$div(
#'       id = "sigpolytope_splash",
#'       tags$div(
#'         class = "card",
#'         tags$h2("SigPolytope — Geometric Omic Signatures"),
#'         tags$p("Welcome to the geometric dashboard for multi-omic regulatory circuitries."),
#'         tags$div(class = "hint", "Click anywhere to start.")
#'       )
#'     ),
#'     
#'     tabItems(
#'         
#'         # ---------------------------------------------------
#'         # 1. OVERVIEW & CONCEPT
#'         # ---------------------------------------------------
#'         tabItem(
#'           tabName = "overview",
#'           fluidRow(
#'             box(
#'               width = 12,
#'               title = "What this dashboard represents",
#'               status = "primary",
#'               solidHeader = TRUE,
#'               collapsible = TRUE,
#'               
#'               p(
#'                 "SigPolytope is an interactive environment for exploring multi-omic regulatory circuitries ",
#'                 "as geometric objects. Instead of treating signatures as isolated variables or gene lists, ",
#'                 "each circuitry is modeled as a structured multidimensional entity whose internal organization ",
#'                 "is preserved and visualized."
#'               ),
#'               
#'               p("Within this dashboard, users can:"),
#'               
#'               tags$ul(
#'                 tags$li("Explore regulatory circuitries embedded in a latent multidimensional space integrating omic, phenotypic, survival, microenvironmental, and immune dimensions."),
#'                 tags$li("Inspect geometric discordance between signature and interaction components using barycenter distances."),
#'                 tags$li("Assess latent regulatory complexity through convex-hull volumes and asymmetry."),
#'                 tags$li("Navigate metabolic superfamilies, pathways, and cell-death programs across geometric regimes."),
#'                 tags$li("Interactively visualize each circuitry as a paired 3D polytope representation.")
#'               ),
#'               
#'               p(
#'                 "Together, these representations provide an intuitive geometric perspective on how ",
#'                 "metabolic regulatory circuitries are organized, coordinated, or divergent across cancers."
#'               )
#'             )
#'           ),
#'           
#'           fluidRow(
#'             box(
#'               width = 6,
#'               title = "Quick start",
#'               status = "info",
#'               solidHeader = TRUE,
#'               tags$ol(
#'                 tags$li("Open the 'Regulatory Circuitries (3D)' tab to access the geometric explorer."),
#'                 tags$li("Apply filters to restrict circuitries by cancer type and metabolic context."),
#'                 tags$li("Select a Circuitries_id to visualize its dual polytope representation."),
#'                 tags$li("Use the summary panel to interpret concordance, discordance, and complexity."),
#'                 tags$li("Switch to the network tab to inspect regulatory interactions in graph form.")
#'               )
#'             ),
#'             
#'             box(
#'               width = 6,
#'               title = "Target-centric access",
#'               status = "info",
#'               solidHeader = TRUE,
#'               p(
#'                 "If your analysis starts from a specific gene, miRNA, lncRNA, or transcript isoform, ",
#'                 "use the OncoMetabolismGPS target-search tools to identify relevant signatures, ",
#'                 "then transition to SigPolytope to examine their geometric regulatory organization."
#'               )
#'             )
#'           )
#'         ),
#'       
#'       
#'       # ---------------------------------------------------
#'       # 2. REGULATORY CIRCUITRIES (3D) — módulo unificado
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "geometry",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             status = "primary",
#'             solidHeader = TRUE,
#'             title = "3D Geometric Explorer of Regulatory Circuitries",
#'             
#'             # ---- TODO O CONTEÚDO DA ABA ESTÁ AQUI ----
#'             mod_circuitry_polytope_ui("poly3d")
#'           )
#'         )
#'       ),
#'       
#'       
#'       # ---------------------------------------------------
#'       # 3. NETWORK VIEWER (já existente)
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "reg_circ",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "Regulatory circuitries — network and textual summary",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             mod_meaningful_interaction_ui("Mean_int")
#'           )
#'         )
#'       ),
#'       
#'       # ---------------------------------------------------
#'       # 4. ABOUT & CITATION
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "about",
#'         
#'         # Bloco 1: How to cite (em cima)
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "How to cite this geometric framework",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             
#'             h3("How to Cite"),
#'             p("If you use this geometric representation or any output generated by this dashboard, please cite:"),
#'             
#'             div(
#'               style = "margin-left:10px;",
#'               p(
#'                 "Nogueira, H. A. C.; Souza, E. R.; Lopes, V. S.; Medina-Acosta, E. (2025). ",
#'                 em("A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer."),
#'                 " Preprint. DOI: ",
#'                 a(
#'                   "10.1101/2025.11.15.688631",
#'                   href   = "https://doi.org/10.1101/2025.11.15.688631",
#'                   target = "_blank"
#'                 )
#'               )
#'             ),
#'             
#'             h4("BibTeX"),
#'             tags$pre(
#'               "@article{Nogueira2025MultiOmicAtlas,
#'               title={A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer},
#'               author={Nogueira, Higor A. C. and Souza, E. R. and Lopes, V. S. and Medina-Acosta, E.},
#'               year={2025},
#'               journal={Preprint},
#'               doi={10.1101/2025.11.15.688631}
#'               }"
#'             )
#'           )
#'         ),
#'         
#'         # Bloco 2: About us (embaixo)
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "About us",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             mod_developers_ui("Developers")
#'           )
#'         )
#'       ),
#'       
#'       
#'       
#'       
#'       
#'       # ---------------------------------------------------
#'       # 5. GITHUB LINK
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "github_link",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "GitHub repository",
#'             status = "info",
#'             solidHeader = TRUE,
#'             p("Click the icon below to open the project repository:"),
#'             tags$a(
#'               href="https://github.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS",
#'               target="_blank",
#'               icon("github"), " BioCancerInformatics / Multi-omic-Oncometabolism-GPS"
#'             )
#'           )
#'         )
#'       )
#'     ),
#'     
#'     # ---- JS para abrir GitHub em nova aba ----
#'     tags$script(HTML("
#'       $(document).on('click', 'a[data-toggle=\"tab\"][data-value=\"github_link\"]', function(e){
#'         e.preventDefault();
#'         window.open('https://github.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS', '_blank');
#'       });
#'     "))
#'   )
#' )
#' 


# =========================================================
# UI - Geometric Multidimensional Representation Dashboard
# =========================================================

# ---- Head tags (CSS + JS) definidos fora do ui ----
APP_HEAD_TAGS <- shiny::singleton(
  tags$head(
    
    # Título da aba do navegador
    tags$title("SigPolytope Shiny"),
    
    tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f7f7f7; }
      .box { border-radius: 10px; }
      .small-box { border-radius: 10px !important; }

      /* =========================
         Loading overlay
         ========================= */
      #app-loading-overlay {
        position: fixed;
        inset: 0;
        background: rgba(255,255,255,0.97);
        z-index: 99999;
        display: flex;
        align-items: center;
        justify-content: center;
        flex-direction: column;
        font-family: Arial, sans-serif;
      }

      #app-loading-title {
        font-size: 26px;
        font-weight: 700;
        color: #111;
        margin-bottom: 10px;
      }

      #app-loading-subtitle {
        font-size: 14px;
        color: #555;
        max-width: 720px;
        text-align: center;
        line-height: 1.5;
        margin-bottom: 18px;
      }

      .spinner {
        width: 44px;
        height: 44px;
        border: 5px solid #ddd;
        border-top: 5px solid #111;
        border-radius: 50%;
        animation: spin 0.9s linear infinite;
      }

      @keyframes spin {
        to { transform: rotate(360deg); }
      }

      /* =========================
         Splash overlay
         ========================= */
      #sigpolytope_splash {
        position: fixed;
        inset: 0;
        z-index: 99998;
        background: rgba(0,0,0,0.70);
        display: none;
        align-items: center;
        justify-content: center;
        padding: 24px;
      }

      #sigpolytope_splash .card {
        max-width: 980px;
        width: 100%;
        background: #ffffff;
        border-radius: 14px;
        padding: 22px 24px;
        box-shadow: 0 10px 30px rgba(0,0,0,0.25);
      }

      #sigpolytope_splash h2 {
        margin-top: 0;
        font-weight: 700;
      }

      #sigpolytope_splash p {
        font-size: 15px;
        line-height: 1.45;
        color: #222;
      }

      #sigpolytope_splash .hint {
        margin-top: 14px;
        font-size: 13px;
        color: #666;
      }
    ")),
    
    tags$script(HTML("
      (function(){
        function hideSplash(){
          var el = document.getElementById('sigpolytope_splash');
          if(el){ el.style.display = 'none'; }
          document.removeEventListener('click', hideSplash, true);
          document.removeEventListener('keydown', hideSplash, true);
          document.removeEventListener('wheel', hideSplash, true);
          document.removeEventListener('touchstart', hideSplash, true);
        }

        function showSplash(){
          var el = document.getElementById('sigpolytope_splash');
          if(el){ el.style.display = 'flex'; }
          document.addEventListener('click', hideSplash, true);
          document.addEventListener('keydown', hideSplash, true);
          document.addEventListener('wheel', hideSplash, true);
          document.addEventListener('touchstart', hideSplash, true);
        }

        document.addEventListener('DOMContentLoaded', function(){
          setTimeout(function(){
            var load = document.getElementById('app-loading-overlay');
            if(load){ load.style.display = 'none'; }
            showSplash();
          }, 5000);
        });
      })();
    "))
  )
)

# =========================================================
# UI principal
# =========================================================

ui <- dashboardPage(
  
  title = "SigPolytope Shiny – OncoMetabolismGPS",
  skin  = "black",
  
  dashboardHeader(
    title = tagList(
      span("SigPolytope Shiny", style = "font-weight:600;"),
      HTML("<span style='font-size:11px; margin-left:4px;'>(OncoMetabolismGPS – Geometry)</span>")
    ),
    titleWidth = 260
  ),
  
  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_menu",
      
      menuItem("Overview & Concept", tabName = "overview", icon = icon("home")),
      menuItem("Regulatory Circuitries (3D)", tabName = "geometry", icon = icon("cube")),
      menuItem("Regulatory Circuitries (network)", tabName = "reg_circ", icon = icon("share-alt")),
      menuItem("About & Citation", tabName = "about", icon = icon("info-circle")),
      
      menuItem(
        "GitHub",
        icon = icon("github"),
        href = "https://github.com/BioCancerInformatics/SigPolytope/tree/main",
        newtab = TRUE
      )
    )
  ),
  
  dashboardBody(
    
    APP_HEAD_TAGS,
    
    # Loading overlay
    tags$div(
      id = "app-loading-overlay",
      tags$div(id = "app-loading-title", "SigPolytope Shiny"),
      tags$div(
        id = "app-loading-subtitle",
        HTML("A geometric atlas of multi-omic regulatory circuitries.")
      ),
      tags$div(class = "spinner"),
      tags$div(style="margin-top:12px;color:#777;font-size:12px;", "Loading…")
    ),
    
    # Splash screen
    tags$div(
      id = "sigpolytope_splash",
      tags$div(
        class = "card",
        tags$h2("SigPolytope — Geometric Omic Signatures"),
        tags$p("Welcome to the geometric dashboard for multi-omic regulatory circuitries."),
        tags$div(class = "hint", "Click anywhere to start.")
      )
    ),
    
    tabItems(
      
      # ---------------------------------------------------
      # OVERVIEW
      # ---------------------------------------------------
      tabItem(
        tabName = "overview",
        fluidRow(
          box(
            width = 12,
            title = "What this dashboard represents",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            
            p("SigPolytope is an interactive environment for exploring multi-omic regulatory circuitries as geometric objects."),
            p("Each circuitry is modeled as a structured multidimensional entity whose internal organization is preserved and visualized."),
            
            tags$ul(
              tags$li("Latent multidimensional embeddings"),
              tags$li("Geometric discordance and asymmetry"),
              tags$li("Convex-hull volume complexity"),
              tags$li("Metabolic and immune coordination"),
              tags$li("Paired 3D polytope visualization")
            )
          )
        )
      ),
      
      # ---------------------------------------------------
      # 3D GEOMETRY
      # ---------------------------------------------------
      tabItem(
        tabName = "geometry",
        fluidRow(
          box(
            width = 12,
            title = "3D Geometric Explorer of Regulatory Circuitries",
            status = "primary",
            solidHeader = TRUE,
            
            # Módulo carregado como no app antigo
            mod_circuitry_polytope_ui("poly3d")
          )
        )
      ),
      
      # ---------------------------------------------------
      # NETWORK
      # ---------------------------------------------------
      tabItem(
        tabName = "reg_circ",
        fluidRow(
          box(
            width = 12,
            title = "Regulatory circuitries — network and textual summary",
            status = "primary",
            solidHeader = TRUE,
            
            mod_meaningful_interaction_ui("Mean_int")
          )
        )
      ),
      
      # ---------------------------------------------------
      # ABOUT
      # ---------------------------------------------------
      tabItem(
        tabName = "about",
        
        fluidRow(
          box(
            width = 12,
            title = "How to cite this geometric framework",
            status = "primary",
            solidHeader = TRUE,
            
            p("If you use this framework, please cite:"),
            tags$pre(
              "@article{Nogueira2025MultiOmicAtlas,
title={A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer},
author={Nogueira, Higor A. C. and Souza, E. R. and Lopes, V. S. and Medina-Acosta, E.},
year={2025},
doi={10.1101/2025.11.15.688631}
}"
            )
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "About us",
            status = "primary",
            solidHeader = TRUE,
            
            mod_developers_ui("Developers")
          )
        )
      )
    )
  )
)







#' 
#' # =========================================================
#' # UI - Geometric Multidimensional Representation Dashboard
#' # =========================================================
#' 
#' # ---- Head tags (CSS + JS) definidos fora do ui para reduzir custo por sessão ----
#' APP_HEAD_TAGS <- shiny::singleton(
#'   tags$head(
#'     
#'     # ✅ Corrige o título do navegador (Chrome tab)
#'     tags$title("SigPolytope Shiny"),
#'     
#'     tags$style(HTML("
#'       .content-wrapper, .right-side { background-color: #f7f7f7; }
#'       .box { border-radius: 10px; }
#'       .small-box { border-radius: 10px !important; }
#' 
#'       /* =========================
#'          Loading overlay (5s)
#'          ========================= */
#'       #app-loading-overlay {
#'         position: fixed; top: 0; left: 0;
#'         width: 100vw; height: 100vh;
#'         background: rgba(255,255,255,0.97);
#'         z-index: 99999;
#'         display: flex;
#'         align-items: center;
#'         justify-content: center;
#'         flex-direction: column;
#'         font-family: Arial, sans-serif;
#'       }
#' 
#'       #app-loading-title {
#'         font-size: 26px;
#'         font-weight: 700;
#'         color: #111;
#'         margin-bottom: 10px;
#'       }
#' 
#'       #app-loading-subtitle {
#'         font-size: 14px;
#'         color: #555;
#'         max-width: 720px;
#'         text-align: center;
#'         line-height: 1.5;
#'         margin-bottom: 18px;
#'       }
#' 
#'       .spinner {
#'         width: 44px; height: 44px;
#'         border: 5px solid #ddd;
#'         border-top: 5px solid #111;
#'         border-radius: 50%;
#'         animation: spin 0.9s linear infinite;
#'       }
#' 
#'       @keyframes spin {
#'         0% { transform: rotate(0deg); }
#'         100% { transform: rotate(360deg); }
#'       }
#' 
#'       /* =========================
#'          Splash overlay (click-to-dismiss)
#'          ========================= */
#'       #sigpolytope_splash {
#'         position: fixed;
#'         inset: 0;
#'         z-index: 99998; /* abaixo do loading (99999) */
#'         background: rgba(0,0,0,0.70);
#'         display: none;  /* aparece após loading */
#'         align-items: center;
#'         justify-content: center;
#'         padding: 24px;
#'       }
#' 
#'       #sigpolytope_splash .card {
#'         max-width: 980px;
#'         width: 100%;
#'         background: #ffffff;
#'         border-radius: 14px;
#'         padding: 22px 24px;
#'         box-shadow: 0 10px 30px rgba(0,0,0,0.25);
#'       }
#' 
#'       #sigpolytope_splash h2 {
#'         margin-top: 0;
#'         font-weight: 700;
#'       }
#' 
#'       #sigpolytope_splash p {
#'         font-size: 15px;
#'         line-height: 1.45;
#'         margin-bottom: 10px;
#'         color: #222;
#'       }
#' 
#'       #sigpolytope_splash .hint {
#'         margin-top: 14px;
#'         font-size: 13px;
#'         color: #666;
#'       }
#'     ")),
#'     
#'     # (mantive o seu JS original; se quiser eu também posso te dar a versão mais leve)
#'     tags$script(HTML("
#'       (function(){
#' 
#'         function hideSplash(){
#'           var el = document.getElementById('sigpolytope_splash');
#'           if(el){ el.style.display = 'none'; }
#'           document.removeEventListener('click', hideSplash, true);
#'           document.removeEventListener('keydown', hideSplash, true);
#'           document.removeEventListener('wheel', hideSplash, true);
#'           document.removeEventListener('touchstart', hideSplash, true);
#'         }
#' 
#'         function showSplash(){
#'           var el = document.getElementById('sigpolytope_splash');
#'           if(el){ el.style.display = 'flex'; }
#'           // some no primeiro clique/interação
#'           document.addEventListener('click', hideSplash, true);
#'           document.addEventListener('keydown', hideSplash, true);
#'           document.addEventListener('wheel', hideSplash, true);
#'           document.addEventListener('touchstart', hideSplash, true);
#'         }
#' 
#'         function hideLoadingAndShowSplash(){
#'           var load = document.getElementById('app-loading-overlay');
#'           if(load){ load.style.display = 'none'; }
#'           showSplash();
#'         }
#' 
#'         // mantém loading por 5s e depois mostra splash
#'         document.addEventListener('DOMContentLoaded', function(){
#'           setTimeout(hideLoadingAndShowSplash, 5000);
#'         });
#' 
#'       })();
#'     "))
#'   )
#' )
#' 
#' ui <- dashboardPage(
#'   # ✅ Também define o título aqui (ajuda e não atrapalha)
#'   title = "SigPolytope Shiny – OncoMetabolismGPS",
#'   
#'   skin = "black",
#'   
#'   dashboardHeader(
#'     title = tagList(
#'       span("SigPolytope Shiny", style = "font-weight:600;"),
#'       HTML('<span style="font-size:11px; margin-left:4px;">(OncoMetabolismGPS – Geometry)</span>')
#'     ),
#'     titleWidth = 260
#'   ),
#'   
#'   dashboardSidebar(
#'     width = 260,
#'     sidebarMenu(
#'       id = "main_menu",
#'       
#'       menuItem("Overview & Concept", tabName = "overview", icon = icon("home")),
#'       menuItem("Regulatory Circuitries (3D)", tabName = "geometry", icon = icon("cube")),
#'       menuItem("Regulatory Circuitries (network)", tabName = "reg_circ", icon = icon("share-alt")),
#'       menuItem("About & Citation", tabName = "about", icon = icon("info-circle")),
#'       
#'       # ---- link direto (remove aba + remove JS jQuery) ----
#'       menuItem(
#'         "GitHub",
#'         icon = icon("github"),
#'         href = "https://github.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS",
#'         newtab = TRUE
#'       )
#'     )
#'   ),
#'   
#'   dashboardBody(
#'     
#'     # ------------------------------------------------------------------
#'     # CSS + JS (singleton, pré-definidos)
#'     # ------------------------------------------------------------------
#'     APP_HEAD_TAGS,
#'     
#'     # ------------------------------------------------------------------
#'     # Loading overlay (aparece imediatamente)
#'     # ------------------------------------------------------------------
#'     tags$div(
#'       id = "app-loading-overlay",
#'       tags$div(id = "app-loading-title", "SigPolytope Shiny"),
#'       tags$div(
#'         id = "app-loading-subtitle",
#'         HTML("A geometric atlas of multi-omic regulatory circuitries.")
#'       ),
#'       tags$div(class = "spinner"),
#'       tags$div(style="margin-top:12px;color:#777;font-size:12px;", "Loading…")
#'     ),
#'     
#'     # ------------------------------------------------------------------
#'     # Splash / mensagem inicial (aparece após 5s; some no clique)
#'     # ------------------------------------------------------------------
#'     tags$div(
#'       id = "sigpolytope_splash",
#'       tags$div(
#'         class = "card",
#'         tags$h2("SigPolytope — Geometric Omic Signatures"),
#'         tags$p("Welcome to the geometric dashboard for multi-omic regulatory circuitries."),
#'         tags$div(class = "hint", "Click anywhere to start.")
#'       )
#'     ),
#'     
#'     tabItems(
#'       
#'       # ---------------------------------------------------
#'       # 1. OVERVIEW & CONCEPT
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "overview",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "What this dashboard represents",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             collapsible = TRUE,
#'             
#'             p(
#'               "SigPolytope is an interactive environment for exploring multi-omic regulatory circuitries ",
#'               "as geometric objects. Instead of treating signatures as isolated variables or gene lists, ",
#'               "each circuitry is modeled as a structured multidimensional entity whose internal organization ",
#'               "is preserved and visualized."
#'             ),
#'             
#'             p("Within this dashboard, users can:"),
#'             
#'             tags$ul(
#'               tags$li("Explore regulatory circuitries embedded in a latent multidimensional space integrating omic, phenotypic, survival, microenvironmental, and immune dimensions."),
#'               tags$li("Inspect geometric discordance between signature and interaction components using barycenter distances."),
#'               tags$li("Assess latent regulatory complexity through convex-hull volumes and asymmetry."),
#'               tags$li("Navigate metabolic superfamilies, pathways, and cell-death programs across geometric regimes."),
#'               tags$li("Interactively visualize each circuitry as a paired 3D polytope representation.")
#'             ),
#'             
#'             p(
#'               "Together, these representations provide an intuitive geometric perspective on how ",
#'               "metabolic regulatory circuitries are organized, coordinated, or divergent across cancers."
#'             )
#'           )
#'         ),
#'         
#'         fluidRow(
#'           box(
#'             width = 6,
#'             title = "Quick start",
#'             status = "info",
#'             solidHeader = TRUE,
#'             tags$ol(
#'               tags$li("Open the 'Regulatory Circuitries (3D)' tab to access the geometric explorer."),
#'               tags$li("Apply filters to restrict circuitries by cancer type and metabolic context."),
#'               tags$li("Select a Circuitries_id to visualize its dual polytope representation."),
#'               tags$li("Use the summary panel to interpret concordance, discordance, and complexity."),
#'               tags$li("Switch to the network tab to inspect regulatory interactions in graph form.")
#'             )
#'           ),
#'           
#'           box(
#'             width = 6,
#'             title = "Target-centric access",
#'             status = "info",
#'             solidHeader = TRUE,
#'             p(
#'               "If your analysis starts from a specific gene, miRNA, lncRNA, or transcript isoform, ",
#'               "use the OncoMetabolismGPS target-search tools to identify relevant signatures, ",
#'               "then transition to SigPolytope to examine their geometric regulatory organization."
#'             )
#'           )
#'         )
#'       ),
#'       
#'       # ---------------------------------------------------
#'       # 2. REGULATORY CIRCUITRIES (3D) — sem lazy-load
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "geometry",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             status = "primary",
#'             solidHeader = TRUE,
#'             title = "3D Geometric Explorer of Regulatory Circuitries",
#'             
#'             # ✅ SEM conditionalPanel: módulo existe desde o início (como no Shiny antigo)
#'             mod_circuitry_polytope_ui("poly3d")
#'           )
#'         )
#'       ),
#'       
#'       # ---------------------------------------------------
#'       # 3. NETWORK VIEWER — sem lazy-load
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "reg_circ",
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "Regulatory circuitries — network and textual summary",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             
#'             # ✅ SEM conditionalPanel: módulo existe desde o início (como no Shiny antigo)
#'             mod_meaningful_interaction_ui("Mean_int")
#'           )
#'         )
#'       ),
#'       
#'       # ---------------------------------------------------
#'       # 4. ABOUT & CITATION — sem lazy-load (mantém igual)
#'       # ---------------------------------------------------
#'       tabItem(
#'         tabName = "about",
#'         
#'         # Bloco 1: How to cite (em cima)
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "How to cite this geometric framework",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             
#'             h3("How to Cite"),
#'             p("If you use this geometric representation or any output generated by this dashboard, please cite:"),
#'             
#'             div(
#'               style = "margin-left:10px;",
#'               p(
#'                 "Nogueira, H. A. C.; Souza, E. R.; Lopes, V. S.; Medina-Acosta, E. (2025). ",
#'                 em("A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer."),
#'                 " Preprint. DOI: ",
#'                 a(
#'                   "10.1101/2025.11.15.688631",
#'                   href   = "https://doi.org/10.1101/2025.11.15.688631",
#'                   target = "_blank"
#'                 )
#'               )
#'             ),
#'             
#'             h4("BibTeX"),
#'             tags$pre(
#'               "@article{Nogueira2025MultiOmicAtlas,
#'               title={A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer},
#'               author={Nogueira, Higor A. C. and Souza, E. R. and Lopes, V. S. and Medina-Acosta, E.},
#'               year={2025},
#'               journal={Preprint},
#'               doi={10.1101/2025.11.15.688631}
#'               }"
#'             )
#'           )
#'         ),
#'         
#'         # Bloco 2: About us (embaixo)
#'         fluidRow(
#'           box(
#'             width = 12,
#'             title = "About us",
#'             status = "primary",
#'             solidHeader = TRUE,
#'             mod_developers_ui("Developers")
#'           )
#'         )
#'       )
#'     )
#'   )
#' )

