# Carregar os arquivos UI e Server
source("ui.R")
source("server.R")

# Iniciar a aplicação
shinyApp(ui, server)

# Run the app
shiny::runApp()


