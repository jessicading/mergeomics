library(shiny)
setwd("/srv/shiny-server/Network_drug_repositioning")
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Network_drug_repositioning")
source("helpers.r")

ui <- fluidPage(
  headerPanel("PharmOmics drug repositioning system"),
  a("PharmOmics Homepage", href="http://mergeomics.research.idre.ucla.edu/PharmOmics/"),
  
  sidebarPanel(
    
    #textInput("mailaddress", "your email address:", ""),
    selectizeInput(
      'system', 'Select network system you want to test repositioning',
      choices = c("Human_liver_network","Rat_liver_network")
    ),
    textAreaInput("Genes", "Input target gene you want to test repositioning, separated by fixed deliminators. Hyperlipidemia from MergeOmics pipeline is used as example here", value = "FASN\nFDFT1\nSQLE\nSC4MOL\nINSIG1\nLSS\nDHCR7\nMGLL\nACLY\nIDI1\nHMGCR\nMVD\nPNPLA5\nPMVK\nNSDHL\nMUM1\nHMGCS1\nDNTT\nELOVL6\nTMEM97\nPEX11A\nFDPS\nDHCR24\nNFE2\nFADS2\nVNN1\nELOVL5\nALDH3A2\nTKT\nPKLR\nPLTP\nGPAM\nPSTPIP2\nMGST3\nHSD17B7\nRDH11\nPAOX\nCOL6A3\nACSS2\nSCD1\nERMP1\nCIDEC\nFRMD4B\nELOVL5\nMOGAT1\nIRF5\nMID1IP1\nTHRSP\nDUSP9\nSTARD4\nSLC2A4\nCOL1A1\nTM7SF2\nACACA\nFERMT3\nDCI\nGNGT1\nAPOH\nACAA1B", width = NULL, height = NULL,
                  cols = 20, rows = 6, placeholder = NULL, resize = NULL),
    actionButton("goButton",label = "Submit job"),
    uiOutput("download_result")
  ), mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel("Progress update",verbatimTextOutput("message")),
      tabPanel("Network repositioning result", tableOutput("analysisresult"))
      #tabPanel("Jaccard and z-score correlation plot",plotOutput("correlation_result"))
    )
  )
)

server <- function(input, output) {
  analysisoutput <- reactiveValues()
  analysisoutput$tableresult <- NULL
  analysisoutput$message <- "Ready for analysis, it will take about 5 minutes to complete analysis"
  
  observe({
    # Take a dependency on input$goButton
    if (input$goButton == 0)
      return(NULL)
    # Use isolate() to avoid dependency on input$goButton
    isolate({
      Genes <- unlist(strsplit(input$Genes,"\n|\t|,| "))
      withProgress({
      setProgress(message = "Checking genes")
      if(input$system %in% "Human_liver_network"){
        diseasegenes <- Genes[Genes %in% V(net_human)$name]
        if(length(diseasegenes) <= 5){
          Genes <- RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes]
          diseasegenes <- Genes[Genes %in% V(net_human)$name]
        }
        if(length(diseasegenes) <= 5){
          setProgress(message = "please check your input format or your genes are not on the network, at least 5 genes are required")
          analysisoutput$message <- "Error, please check your input format or your genes are not on the network, at least 5 genes are required"
        }else{
          setProgress(message = "Start analysis, will take minutes")
          net <- net_human
          DOBtable <- DOBtable_human
          alldistancetable <- alldistancetable_human
          allchangegenelist <- allchangegenelist_human
          DEGlist_up <- human_DEGlist_up
          DEGlist_down <- human_DEGlist_down
          
          Diseasedeg <- degree(net, v = diseasegenes, mode = "all",
                               loops = TRUE, normalized = FALSE)
          
          DiseaseDEGtable <- as.data.frame(table(Diseasedeg), stringsAsFactors = FALSE)
          DiseaseDEGtable$Diseasedeg <- as.numeric(DiseaseDEGtable$Diseasedeg)
          DiseaseDEGtable$Freq <- as.numeric(DiseaseDEGtable$Freq)
          
          
          for(i in 1:nrow(DiseaseDEGtable)){
            genenumber <- DiseaseDEGtable$Freq[i]
            DEGnumber <- DiseaseDEGtable$Diseasedeg[i]
            
            if(DEGnumber >= 31 & DEGnumber < 33){DEGnumber <- 31:32
            }else if(DEGnumber >= 33 & DEGnumber < 36){DEGnumber <- 33:35
            }else if(DEGnumber >= 36 & DEGnumber < 40){DEGnumber <- 36:39
            }else if(DEGnumber >= 40 & DEGnumber < 45){DEGnumber <- 40:44
            }else if(DEGnumber >= 45 & DEGnumber < 51){DEGnumber <- 45:50
            }else if(DEGnumber >= 51){DEGnumber <- 51:201}
            genepools <- DOBtable$gene[DOBtable$degree %in% DEGnumber]
            sampledgenes <- replicate(1000,sample(genepools,genenumber,replace = TRUE))
            if(!exists("finalsample")){
              finalsample <- sampledgenes
            }else{
              finalsample <- rbind(finalsample,sampledgenes)
            }
          }
          if(class(finalsample) %in% "character"){
            randomgenes_disease <- as.list(finalsample)
          }else{
            randomgenes_disease <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
            
          }
          resd <- NULL
          Jaccard_resd <- NULL
          alldrugs <- names(allchangegenelist)
          for(j in 1:length(alldrugs )){
            changegenes <- allchangegenelist[[alldrugs[j]]]
            Drugdeg <- degree(net, v = changegenes, mode = "all",
                              loops = TRUE, normalized = FALSE)
            randomgenes_drug <- list()
            DrugDEGtable <- as.data.frame(table(Drugdeg), stringsAsFactors = FALSE)
            DrugDEGtable[,2] <- as.numeric(DrugDEGtable[,2])
            DrugDEGtable[,1] <- as.numeric(DrugDEGtable[,1])
            
            rm(finalsample)
            for(i in 1:nrow(DrugDEGtable)){
              genenumber <- DrugDEGtable$Freq[i]
              DEGnumber <- DrugDEGtable$Drugdeg[i]
              
              if(DEGnumber >= 31 & DEGnumber < 33){DEGnumber <- 31:32
              }else if(DEGnumber >= 33 & DEGnumber < 36){DEGnumber <- 33:35
              }else if(DEGnumber >= 36 & DEGnumber < 40){DEGnumber <- 36:39
              }else if(DEGnumber >= 40 & DEGnumber < 45){DEGnumber <- 40:44
              }else if(DEGnumber >= 45 & DEGnumber < 51){DEGnumber <- 45:50
              }else if(DEGnumber >= 51){DEGnumber <- 51:201}
              genepools <- DOBtable$gene[DOBtable$degree %in% DEGnumber]
              sampledgenes <- replicate(1000,sample(genepools,genenumber,replace = TRUE))
              if(!exists("finalsample")){
                finalsample <- sampledgenes
              }else{
                finalsample <- rbind(finalsample,sampledgenes)
              }
            }
            if(class(finalsample) %in% "character"){
              randomgenes_drug <- as.list(finalsample)
            }else{
              randomgenes_drug <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
              
            }

            
            allscores <- NULL
            for(k in 1:1000){
              
              if(length(changegenes) == 1){
                allscores[k] <- mean(alldistancetable[randomgenes_drug[[k]],randomgenes_disease[[k]]], na.rm = T)
              }else{
                allscores[k] <- mean(apply( alldistancetable[randomgenes_drug[[k]],randomgenes_disease[[k]]], 1, min, na.rm = T ), na.rm = T)
              }
            }
            if(length(changegenes) == 1){
              res <- mean(alldistancetable[changegenes,diseasegenes], na.rm = T)
            }else{
              res <- mean(apply( alldistancetable[changegenes,diseasegenes], 1,min, na.rm = T ), na.rm = T)
            }
            resd[alldrugs[j]] <- (res-mean(allscores, na.rm = T))/sd(allscores,na.rm = T)
            Jaccard_resd[alldrugs[j]] <- Jaccard2(unique(c(DEGlist_up[[alldrugs[j]]], DEGlist_down[[alldrugs[j]]])),Genes)
          }
          analysisoutput$message <-  setProgress(message = "You may check score table and download result")
          analysisoutput$tableresult <- data.frame(drugname = names(resd),z_score = resd, Jaccard_score = Jaccard_resd)
        }
      }else{
        diseasegenes <- Genes[Genes %in% V(net_mouse)$name]
        if(length(diseasegenes) <= 5){
          Genes <- RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes]
          diseasegenes <- Genes[Genes %in% V(net_mouse)$name]
        }
        if(length(diseasegenes) <= 5){
          setProgress(message = "please check your input format or your genes are not on the network, at least 5 genes are required")
          analysisoutput$message <- "Error, please check your input format or your genes are not on the network, at least 5 genes are required"
        }else{
          setProgress(message = "Start analysis, will take minutes")
          net <- net_mouse
          DOBtable <- DOBtable_mouse
          alldistancetable <- alldistancetable_mouse
          allchangegenelist <- allchangegenelist_mouse
          DEGlist_up <- Rat_DEGlist_up
          DEGlist_down <- Rat_DEGlist_down
          
          Diseasedeg <- degree(net, v = diseasegenes, mode = "all",
                               loops = TRUE, normalized = FALSE)
          
          DiseaseDEGtable <- as.data.frame(table(Diseasedeg), stringsAsFactors = FALSE)
          DiseaseDEGtable$Diseasedeg <- as.numeric(DiseaseDEGtable$Diseasedeg)
          DiseaseDEGtable$Freq <- as.numeric(DiseaseDEGtable$Freq)
          
          rm(finalsample)
          for(i in 1:nrow(DiseaseDEGtable)){
            genenumber <- DiseaseDEGtable$Freq[i]
            DEGnumber <- DiseaseDEGtable$Diseasedeg[i]
            
            if(DEGnumber >= 31 & DEGnumber < 33){DEGnumber <- 31:32
            }else if(DEGnumber >= 33 & DEGnumber < 36){DEGnumber <- 33:35
            }else if(DEGnumber >= 36 & DEGnumber < 40){DEGnumber <- 36:39
            }else if(DEGnumber >= 40 & DEGnumber < 45){DEGnumber <- 40:44
            }else if(DEGnumber >= 45 & DEGnumber < 51){DEGnumber <- 45:50
            }else if(DEGnumber >= 51){DEGnumber <- 51:201}
            genepools <- DOBtable$gene[DOBtable$degree %in% DEGnumber]
            sampledgenes <- replicate(1000,sample(genepools,genenumber,replace = TRUE))
            if(!exists("finalsample")){
              finalsample <- sampledgenes
            }else{
              finalsample <- rbind(finalsample,sampledgenes)
            }
          }
          if(class(finalsample) %in% "character"){
            randomgenes_disease <- as.list(finalsample)
          }else{
            randomgenes_disease <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
            
          }
          
          resd <- NULL
          Jaccard_resd <- NULL
          alldrugs <- names(allchangegenelist)
          for(j in 1:length(alldrugs )){
            changegenes <- allchangegenelist[[alldrugs[j]]]
            Drugdeg <- degree(net, v = changegenes, mode = "all",
                              loops = TRUE, normalized = FALSE)
            randomgenes_drug <- list()
            DrugDEGtable <- as.data.frame(table(Drugdeg), stringsAsFactors = FALSE)
            DrugDEGtable[,2] <- as.numeric(DrugDEGtable[,2])
            DrugDEGtable[,1] <- as.numeric(DrugDEGtable[,1])
            
            rm(finalsample)
            for(i in 1:nrow(DrugDEGtable)){
              genenumber <- DrugDEGtable$Freq[i]
              DEGnumber <- DrugDEGtable$Drugdeg[i]
              
              if(DEGnumber >= 31 & DEGnumber < 33){DEGnumber <- 31:32
              }else if(DEGnumber >= 33 & DEGnumber < 36){DEGnumber <- 33:35
              }else if(DEGnumber >= 36 & DEGnumber < 40){DEGnumber <- 36:39
              }else if(DEGnumber >= 40 & DEGnumber < 45){DEGnumber <- 40:44
              }else if(DEGnumber >= 45 & DEGnumber < 51){DEGnumber <- 45:50
              }else if(DEGnumber >= 51){DEGnumber <- 51:201}
              genepools <- DOBtable$gene[DOBtable$degree %in% DEGnumber]
              sampledgenes <- replicate(1000,sample(genepools,genenumber,replace = TRUE))
              if(!exists("finalsample")){
                finalsample <- sampledgenes
              }else{
                finalsample <- rbind(finalsample,sampledgenes)
              }
            }
            if(class(finalsample) %in% "character"){
              randomgenes_drug <- as.list(finalsample)
            }else{
              randomgenes_drug <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
            }
            
            
            allscores <- NULL
            for(k in 1:1000){
              
              if(length(changegenes) == 1){
                allscores[k] <- mean(alldistancetable[randomgenes_drug[[k]],randomgenes_disease[[k]]], na.rm = T)
              }else{
                allscores[k] <- mean(apply( alldistancetable[randomgenes_drug[[k]],randomgenes_disease[[k]]], 1, min, na.rm = T ), na.rm = T)
              }
            }
            
            if(length(changegenes) == 1){
              res <- mean(alldistancetable[changegenes,diseasegenes], na.rm = T)
            }else{
              res <- mean(apply( alldistancetable[changegenes,diseasegenes], 1,min, na.rm = T ), na.rm = T)
            }
            resd[alldrugs[j]] <- (res-mean(allscores, na.rm = T))/sd(allscores,na.rm = T)
            Jaccard_resd[alldrugs[j]] <- Jaccard2(unique(c(DEGlist_up[[alldrugs[j]]], DEGlist_down[[alldrugs[j]]])),Genes)
          }
          setProgress(message = "You may check score table and download result")
          analysisoutput$tableresult <- data.frame(drugname = names(resd),z_score = resd, Jaccard_score = Jaccard_resd)
        }
      }
      })
      #email <- mime(
      #  To = input$mailaddress,
      #  From = "mergeomics@gmail.com",
      #  Subject = "PharmOmics network repositioning result",
      #  body = "This is an automatic mail, don't reply\n
      #  your drug repositioning result is in attached file with selected network")
      #email <- attach_file(email, file = "/Users/Tsai_Lab/Desktop/Box Sync/LINCS1000_Slicr/lossdrugs.csv", type = "text/plain")
      #email <- attach_file(email, file = "/Users/Tsai_Lab/Desktop/Box Sync/LINCS1000_Slicr/lossdrugs.csv", type = "text/plain")
      
      #send_message(email)
    })
  })
  
  output$analysisresult <- renderTable({
    analysisoutput$tableresult[order(analysisoutput$tableresult$z_score),]
  })
  
  output$correlation_result <- renderPlot({
    if(is.null(analysisoutput$tableresult)){
    }else{
      tableresult <- analysisoutput$tableresult
      corresult <- cor(tableresult$z_score,tableresult$Jaccard_score)
      plot(tableresult$Jaccard_score,tableresult$z_score, las = 1, pch = 20, main = paste0("Jaccard_zscore_correlation, Cor = ",formatC(corresult, format = "e",digits = 2)), xlab = "Jaccard score",ylab = "network_zscore")
    }
  })
  
  output$download_result <- renderUI({
    if(!is.null(analysisoutput$tableresult)) {
      analysisoutput$message <- "Analysis done successfully, you can check results now"
      downloadButton("downloadresult", "Download_result")
    }
  })
  
  output$message <- renderText({analysisoutput$message})
  
  output$downloadresult <- downloadHandler(
    filename = function() {
      "Analysis_result.csv"
    },
    content = function(file) {
      write.csv(analysisoutput$tableresult[order(analysisoutput$tableresult$z_score),], file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui = ui, server = server)
