library(shiny)
Jaccard2 <- function(set1,set2){
  I <- length(intersect(set1,set2))
  return(I/(length(set1)+length(set2)-I))
}
#load("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/All_combined_data_AE_GEO_0.01.rda")
load("/srv/shiny-server/All_combined_data_AE_GEO_0.01.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Shiny_package_new/All_combined_data_AE_GEO_0.01.rda")
ui <- fluidPage(
  headerPanel("PharmOmics drug repositioning system_Jaccard overlapping test"),
  a("PharmOmics Homepage", href="http://mergeomics.research.idre.ucla.edu/PharmOmics/"),
  
  sidebarPanel(
    textAreaInput("Genes_up", "Input upregulated gene you want to test repositioning, separated by fixed deliminators. Hyperlipidemia from MergeOmics pipeline is used as example here", value = "FASN\nFDFT1\nSQLE\nSC4MOL\nINSIG1\nLSS\nDHCR7\nMGLL\nACLY\nIDI1\nHMGCR\nMVD\nPNPLA5\nPMVK\nNSDHL\nMUM1\nHMGCS1\nDNTT\nELOVL6\nTMEM97\nPEX11A\nFDPS\nDHCR24\nNFE2\nFADS2\nVNN1\nELOVL5\nALDH3A2\nTKT\nPKLR\nPLTP\nGPAM\nPSTPIP2\nMGST3\nHSD17B7\nRDH11\nPAOX\nCOL6A3\nACSS2\nSCD1\nERMP1\nCIDEC\nFRMD4B\nELOVL5\nMOGAT1\nIRF5\nMID1IP1\nTHRSP\nDUSP9\nSTARD4\nSLC2A4\nCOL1A1\nTM7SF2\nACACA\nFERMT3\nDCI\nGNGT1\nAPOH\nACAA1B", width = NULL, height = NULL,
                  cols = 20, rows = 6, placeholder = NULL, resize = NULL),
    textAreaInput("Genes_down", "Input downregulated gene you want to test repositioning (optional), separated by fixed deliminators. ", width = NULL, height = NULL,
                  cols = 20, rows = 6, placeholder = NULL, resize = NULL),   
    actionButton("goButton",label = "Submit job"),
    uiOutput("download_result")
  ), mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel("Jaccard repositioning result", tableOutput("analysisresult"))
      #tabPanel("Jaccard and z-score correlation plot",plotOutput("correlation_result"))
    )
  )
)

server <- function(input, output) {
  analysisoutput <- reactiveValues()
  analysisoutput$tableresult <- data.frame(drugname = numeric(0),species = character(0),tissue = character(0), Jaccard_score = numeric(0))

  observe({
    # Take a dependency on input$goButton
    if (input$goButton == 0)
      return(NULL)
    # Use isolate() to avoid dependency on input$goButton
    isolate({
      Genes_up <- unlist(strsplit(input$Genes_up,"\n|\t|,| "))
      Genes_down <- unlist(strsplit(input$Genes_down,"\n|\t|,| "))
      bothind <- ifelse(length(Genes_down) != 0, T,F)
      diseasegenes_up <- Genes_up[Genes_up %in% HUGO_symbols2$`Approved Symbol`]
      if(bothind){diseasegenes_down <- Genes_down[Genes_down %in% HUGO_symbols2$`Approved Symbol`]}
      if(length(diseasegenes_up)/length(diseasegenes_up) < 0.05){
        diseasegenes_up <- unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% diseasegenes_up])
        if(bothind)
        diseasegenes_down <- unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% diseasegenes_down])
      }
      
      drug_all <- NULL
      tissue_all <- NULL
      species_all <- NULL
      Jaccardscore <- NULL
      for(drugs in names(final_list_up)){
        for(species in names(final_list_up[[drugs]])){
          for(tissue in names(final_list_up[[drugs]][[species]])){
            if(tissue %in% "common"){next}
            if(species %in% "Rattus norvigicus"){
              druggenes_up <- unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% final_list_up[[drugs]][[species]][[tissue]]])
              druggenes_down <- unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% final_list_down[[drugs]][[species]][[tissue]]])
            }else if(species %in% "Mus musculus"){
              druggenes_up <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% final_list_up[[drugs]][[species]][[tissue]]])
              druggenes_down <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% final_list_down[[drugs]][[species]][[tissue]]])
            }else{
              druggenes_up <-  final_list_up[[drugs]][[species]][[tissue]]
              druggenes_down <-  final_list_down[[drugs]][[species]][[tissue]]
            }
            drug_all <- c(drug_all,drugs)
            species_all <- c(species_all, species)
            tissue_all <- c(tissue_all, tissue)
            if(bothind){
              Jaccardscore <- c(Jaccardscore, Jaccard2(druggenes_up, diseasegenes_up) + Jaccard2(druggenes_down, diseasegenes_down)-Jaccard2(druggenes_up,diseasegenes_down)-Jaccard2(druggenes_down,diseasegenes_up))
            }else{
              Jaccardscore <- c(Jaccardscore, Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up))
            }
            
          }
        }
      }

      analysisoutput$tableresult <- data.frame(drugname = drug_all,species = species_all,tissue = tissue_all, Jaccard_score = Jaccardscore)
     
    })
  })
  
  output$analysisresult <- renderTable({
    analysisoutput$tableresult[order(analysisoutput$tableresult$Jaccard_score, decreasing = T),]
  })
  
  
  output$download_result <- renderUI({
    if(nrow(analysisoutput$tableresult) != 0) {
      downloadButton("downloadresult", "Download_result")
    }
  })
  
  
  output$downloadresult <- downloadHandler(
    filename = function() {
      "Analysis_result.csv"
    },
    content = function(file) {
      write.csv(analysisoutput$tableresult[order(analysisoutput$tableresult$Jaccard_score, decreasing = T),], file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui = ui, server = server)
