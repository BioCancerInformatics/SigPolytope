mod_developers_ui <- function(id) {
  ns <- NS(id)
  
  # ============================================================
  # CSS INLINE (criado uma única vez, sem sprintf, sem reatividade)
  # ============================================================
  styles <- tags$style(HTML("
.dev-wrap {
  max-width: 960px;
  margin: 0 auto;
  padding: 24px 16px 56px;
  font-family: 'Segoe UI','Roboto',sans-serif;
}
.dev-header-title {
  text-align: center;
  font-size: 1.6rem;
  font-weight: 600;
  color: #1f2933;
  margin-bottom: 4px;
}
.dev-header-subtitle {
  text-align: center;
  color: #6b7280;
  font-size: 0.95rem;
  max-width: 720px;
  margin: 0 auto 24px;
}
.dev-section-title {
  font-size: 1.1rem;
  font-weight: 600;
  color: #1f2933;
  margin: 24px 0 8px;
}
.dev-section-text {
  color: #4b5563;
  font-size: 0.95rem;
  margin-bottom: 12px;
}
.dev-grid {
  display: grid;
  grid-template-columns: repeat(1, minmax(0, 1fr));
  gap: 16px;
}
@media (min-width: 576px) {
  .dev-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
}
@media (min-width: 768px) {
  .dev-grid { grid-template-columns: repeat(4, minmax(0, 1fr)); }
}
.dev-card {
  background: #ffffff;
  border-radius: 12px;
  border: 1px solid #e5e7eb;
  padding: 16px 14px 14px;
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
  gap: 6px;
  box-shadow: 0 4px 10px rgba(15,23,42,0.04);
}
.dev-avatar {
  width: 120px;
  height: 120px;
  border-radius: 999px;
  object-fit: cover;
  border: 2px solid #edf2f7;
  margin-bottom: 6px;
}
.dev-name {
  font-size: 0.98rem;
  font-weight: 600;
  color: #111827;
}
.dev-role {
  font-size: 0.9rem;
  color: #6b7280;
  margin-bottom: 4px;
}
.dev-links {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  gap: 6px;
  margin-top: 4px;
}
.dev-btn {
  display: inline-flex;
  align-items: center;
  gap: 5px;
  padding: 5px 10px;
  border-radius: 999px;
  border: 1px solid #e5e7eb;
  background-color: #f9fafb;
  color: #111827;
  font-size: 0.82rem;
  text-decoration: none;
}
.dev-logo-card {
  background: #ffffff;
  border-radius: 12px;
  border: 1px solid #e5e7eb;
  padding: 18px 16px;
  text-align: center;
  margin-top: 12px;
}
.dev-logo-img {
  width: 140px;
  height: 140px;
  border-radius: 999px;
  object-fit: cover;
  border: 2px solid #edf2f7;
  margin-bottom: 8px;
}
.dev-logo-text {
  font-size: 0.9rem;
  color: #6b7280;
}
.dev-issue-card {
  background: #ffffff;
  border-radius: 12px;
  border: 1px solid #e5e7eb;
  padding: 16px;
  margin-top: 12px;
  font-size: 0.93rem;
  color: #4b5563;
}
.dev-issue-card ul {
  padding-left: 18px;
}
  "))
  
  # ============================================================
  # Dados (estáticos)
  # ============================================================
  people <- list(
    list(
      name = "Higor Almeida Cordeiro Nogueira",
      role = "Researcher",
      img  = "HACN_Photo.jpg",
      rg   = "https://www.researchgate.net/profile/Higor-Cordeiro-Nogueira",
      gh   = "https://github.com/HigorACNogueira",
      li   = "https://linkedin.com/in/higor-almeida-950082255"
    ),
    list(
      name = "Emanuell Rodrigues de Souza",
      role = "Researcher",
      img  = "ESR_Photo.jpg",
      rg   = "https://www.researchgate.net/profile/Emanuell-Rodrigues-De-Souza",
      gh   = "https://github.com/Emanuell-Souza",
      li   = "http://www.linkedin.com/in/emanuell-rodrigues-de-souza-35b40a300"
    ),
    list(
      name = "Victor dos Santos Lopes",
      role = "Researcher",
      img  = "VSL_Photo.png",
      rg   = "https://www.researchgate.net/profile/Victor-Lopes-25",
      li   = "https://www.linkedin.com/in/victor-lopes-880604377"
    ),
    list(
      name = "Enrique Medina-Acosta",
      role = "Team head",
      img  = "EMA_Photo.jpg",
      rg   = "https://www.researchgate.net/profile/Enrique-Medina-Acosta",
      gh   = "https://github.com/quiquemedina",
      loop = "https://loop.frontiersin.org/people/51475/overview"
    )
  )
  
  # ============================================================
  # Card (construção leve, sem tagAppendChild)
  # ============================================================
  person_card <- function(p) {
    
    links <- list(
      if (!is.null(p$rg))
        tags$a(class="dev-btn", href=p$rg, target="_blank", icon("graduation-cap"), "ResearchGate"),
      if (!is.null(p$gh))
        tags$a(class="dev-btn", href=p$gh, target="_blank", icon("github"), "GitHub"),
      if (!is.null(p$li))
        tags$a(class="dev-btn", href=p$li, target="_blank", icon("linkedin"), "LinkedIn"),
      if (!is.null(p$loop))
        tags$a(class="dev-btn", href=p$loop, target="_blank", icon("user-circle"), "Loop")
    )
    
    tags$div(
      class = "dev-card",
      tags$img(
        src = p$img,
        alt = paste("Photo of", p$name),
        class = "dev-avatar",
        loading = "lazy",
        decoding = "async"
      ),
      tags$div(class = "dev-name", p$name),
      tags$div(class = "dev-role", p$role),
      tags$div(class = "dev-links", links)
    )
  }
  
  # ============================================================
  # UI FINAL (100% estático)
  # ============================================================
  tagList(
    styles,
    
    tags$div(
      class = "dev-wrap",
      
      tags$div(class="dev-header-title", "About the project team"),
      tags$p(
        class="dev-header-subtitle",
        "SigPolytope is developed by a research group focused on multi-omic, data-driven approaches to understand cancer metabolism and its clinical implications."
      ),
      
      tags$div(class="dev-section-title", "Team"),
      tags$p(
        class="dev-section-text",
        "These are the people behind the design, curation and implementation of the platform."
      ),
      tags$div(class="dev-grid", lapply(people, person_card)),
      
      tags$div(class="dev-section-title", "Issues and feedback"),
      tags$div(
        class="dev-issue-card",
        tags$p("If you find a problem or have suggestions for improvement, sending a short description already helps a lot."),
        tags$ul(
          tags$li("Where in the app the issue occurred;"),
          tags$li("What you expected to see versus what you observed;"),
          tags$li("A screenshot and error message (if available).")
        ),
        tags$p(
          "You can contact us at ",
          tags$a("higoralmeida1995@gmail.com", href="mailto:higoralmeida1995@gmail.com"),
          "."
        )
      )
    )
  )
}

mod_developers_server <- function(id) {
  moduleServer(id, function(input, output, session) { })
}
