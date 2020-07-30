library(shiny)
Jaccard2 <- function(set1,set2){
  I <- length(intersect(set1,set2))
  return(I/(length(set1)+length(set2)-I))
}

#load("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/PharmOmicsframe.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/All_combined_data_AE_GEO_0.01.rda")#
#load("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/Jaccard_app_database.rda")
#load("/srv/shiny-server/All_combined_data_AE_GEO_0.01.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Shiny_package_new/All_combined_data_AE_GEO_0.01.rda")

#Ratframe <- Combination_everything[Combination_everything$allspecies %in% "Rattus norvegicus",]
#Ratgenesup <- strsplit(Ratframe$allsignatures_up,",")
#Ratgenesdown <- strsplit(Ratframe$allsignatures_down,",")

#Mouseframe <- Combination_everything[Combination_everything$allspecies %in% "Mus musculus",]
#Mousegenesup <- strsplit(Mouseframe$allsignatures_up,",")
#Mousegenesdown <- strsplit(Mouseframe$allsignatures_down,",")

#Humanframe <- Combination_everything[Combination_everything$allspecies %in% "Homo sapiens",]
#Humangenesup <- strsplit(Humanframe$allsignatures_up,",")
#Humangenesdown <- strsplit(Humanframe$allsignatures_down,",")
#save(Ratframe, Ratgenesup,Ratgenesdown,
#     Mouseframe,Mousegenesup,Mousegenesdown,
#     Humanframe,Humangenesup,Humangenesdown, HUGO_symbols2,RAT_symbols2,Mouse_symbols2,
#     file = "/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/Jaccard_app_databasev2.rda")


#load("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/Jaccard_app_databasev2.rda")
load("/srv/shiny-server/Jaccard_repositioning/Jaccard_app_databasev2.rda")
ui <- fluidPage(
  headerPanel("PharmOmics drug repositioning system_Jaccard overlapping test"),
  a("PharmOmics Homepage", href="http://mergeomics.research.idre.ucla.edu/PharmOmics/"),
  
  sidebarPanel(
    textAreaInput("Genes_up", "Input upregulated gene you want to test repositioning, separated by fixed deliminators. Hyperlipidemia from MergeOmics pipeline is used as example here", value = "FASN\nFDFT1\nSQLE\nSC4MOL\nINSIG1\nLSS\nDHCR7\nMGLL\nACLY\nIDI1\nHMGCR\nMVD\nPNPLA5\nPMVK\nNSDHL\nMUM1\nHMGCS1\nDNTT\nELOVL6\nTMEM97\nPEX11A\nFDPS\nDHCR24\nNFE2\nFADS2\nVNN1\nELOVL5\nALDH3A2\nTKT\nPKLR\nPLTP\nGPAM\nPSTPIP2\nMGST3\nHSD17B7\nRDH11\nPAOX\nCOL6A3\nACSS2\nSCD1\nERMP1\nCIDEC\nFRMD4B\nELOVL5\nMOGAT1\nIRF5\nMID1IP1\nTHRSP\nDUSP9\nSTARD4\nSLC2A4\nCOL1A1\nTM7SF2\nACACA\nFERMT3\nDCI\nGNGT1\nAPOH\nACAA1B",
                  width = NULL, height = NULL,
                  cols = 20, rows = 6, placeholder = NULL, resize = NULL),
    textAreaInput("Genes_down", "Input downregulated gene you want to test repositioning (optional), separated by fixed deliminators. ", width = NULL, height = NULL,
                  cols = 20, rows = 6, placeholder = NULL, resize = NULL),   
    actionButton("goButton",label = "Submit job"),
    uiOutput("download_result")
  ), mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel("Jaccard repositioning top rows", tableOutput("analysisresult"))
      #tabPanel("Jaccard and z-score correlation plot",plotOutput("correlation_result"))
    )
  )
)

server <- function(input, output) {
  analysisoutput <- reactiveValues()
  analysisoutput$tableresult <- data.frame("Still running" = character(0))

  observe({
    # Take a dependency on input$goButton
    if (input$goButton == 0)
      return(NULL)
    # Use isolate() to avoid dependency on input$goButton
    isolate({
      Genes_up <- unlist(strsplit(input$Genes_up,"\n|\t|,| "))
      Genes_down <- unlist(strsplit(input$Genes_down,"\n|\t|,| "))
      bothind <- ifelse(length(Genes_down) != 0, T,F)
      
      speciesind <- "Human"
      diseasegenes_up <- Genes_up[Genes_up %in% HUGO_symbols2$`Approved Symbol`]
      if(bothind){diseasegenes_down <- Genes_down[Genes_down %in% HUGO_symbols2$`Approved Symbol`]}
      
      
      if(length(diseasegenes_up)/length(Genes_up) < 0.05){
        species <- "non-Human"
        diseasegenes_up_Rat <- Genes_up[Genes_up %in% RAT_symbols2$rat_symbol]
        if(bothind){diseasegenes_down_Rat <- Genes_down[Genes_down %in% RAT_symbols2$rat_symbol]}
        
        diseasegenes_up <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_up])
        if(bothind){diseasegenes_down <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_down])}
        
        diseasegenes_up_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]
        if(bothind){diseasegenes_down_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]}
        }else{
          diseasegenes_up_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_up])
          if(bothind){diseasegenes_down_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_down])}
          diseasegenes_up_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_up])
          if(bothind){diseasegenes_down_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_down])}
      }

      Jaccardscore <- NULL
      for(r in 1:nrow(Ratframe)){
        druggenes_up <- Ratgenesup[[r]]
        druggenes_down <- Ratgenesdown[[r]]
        if(bothind){
          Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up_Rat) + Jaccard2(druggenes_down, diseasegenes_down_Rat)-Jaccard2(druggenes_up,diseasegenes_down_Rat)-Jaccard2(druggenes_down,diseasegenes_up_Rat)
        }else{
          Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Rat)
        }
      }
      Ratframe$Jaccardscore <- Jaccardscore
      if(bothind){
        rank_within_species <- rep(NA,nrow(Ratframe))
        rank_within_species[which(Ratframe$Jaccardscore > 0)] <- rank( Ratframe$Jaccardscore[Ratframe$Jaccardscore > 0])/nrow( Ratframe[Ratframe$Jaccardscore > 0,])
        rank_within_species[which(Ratframe$Jaccardscore < 0)] <- -rank( -Ratframe$Jaccardscore[Ratframe$Jaccardscore < 0])/nrow( Ratframe[Ratframe$Jaccardscore < 0,])
        Ratframe$rank_within_species <- rank_within_species
      }else{
        Ratframe$rank_within_species <- rank( Ratframe$Jaccardscore)/nrow( Ratframe)
      } 
      
      Jaccardscore <- NULL
      for(r in 1:nrow(Mouseframe)){
        druggenes_up <- Mousegenesup[[r]]
        druggenes_down <- Mousegenesdown[[r]]
        if(bothind){
          Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up_Mouse) + Jaccard2(druggenes_down, diseasegenes_down_Mouse)-Jaccard2(druggenes_up,diseasegenes_down_Mouse)-Jaccard2(druggenes_down,diseasegenes_up_Mouse)
        }else{
          Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Mouse)
        }
      }
      Mouseframe$Jaccardscore <- Jaccardscore
      if(bothind){
        rank_within_species <- rep(NA,nrow(Mouseframe))
        rank_within_species[which(Mouseframe$Jaccardscore > 0)] <- rank( Mouseframe$Jaccardscore[Mouseframe$Jaccardscore > 0])/nrow( Mouseframe[Mouseframe$Jaccardscore > 0,])
        rank_within_species[which(Mouseframe$Jaccardscore < 0)] <- -rank( -Mouseframe$Jaccardscore[Mouseframe$Jaccardscore < 0])/nrow( Mouseframe[Mouseframe$Jaccardscore < 0,])
        Mouseframe$rank_within_species <- rank_within_species
      }else{
        Mouseframe$rank_within_species <- rank( Mouseframe$Jaccardscore)/nrow( Mouseframe)
      }      
      
      Jaccardscore <- NULL
      for(r in 1:nrow(Humanframe)){
        druggenes_up <- Humangenesup[[r]]
        druggenes_down <- Humangenesup[[r]]
        if(bothind){
          Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up) + Jaccard2(druggenes_down, diseasegenes_down)-Jaccard2(druggenes_up,diseasegenes_down)-Jaccard2(druggenes_down,diseasegenes_up)
        }else{
          Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up)
        }
      }
      Humanframe$Jaccardscore <- Jaccardscore
      if(bothind){
        rank_within_species <- rep(NA,nrow(Humanframe))
        rank_within_species[which(Humanframe$Jaccardscore > 0)] <- rank( Humanframe$Jaccardscore[Humanframe$Jaccardscore > 0])/nrow( Humanframe[Humanframe$Jaccardscore > 0,])
        rank_within_species[which(Humanframe$Jaccardscore < 0)] <- -rank( -Humanframe$Jaccardscore[Humanframe$Jaccardscore < 0])/nrow( Humanframe[Humanframe$Jaccardscore < 0,])
        Humanframe$rank_within_species <- rank_within_species
      }else{
        Humanframe$rank_within_species <- rank( Humanframe$Jaccardscore)/nrow( Humanframe)
      }
     
      
      analysisoutput$tableresult <- rbind.data.frame(Humanframe[,-c(6:7)],Mouseframe[,-c(6:7)],Ratframe[,-c(6:7)])
      analysisoutput$tableresult_full <- rbind.data.frame(Humanframe,Mouseframe,Ratframe)
      
    })
  })
  
  output$analysisresult <- renderTable({
    if(nrow(analysisoutput$tableresult) != 0) {
    head(analysisoutput$tableresult[order(analysisoutput$tableresult$Jaccardscore, decreasing = T),], n = 30)
    }
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
      write.csv(analysisoutput$tableresult_full[order(analysisoutput$tableresult_full$Jaccardscore, decreasing = T),], file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui = ui, server = server)
