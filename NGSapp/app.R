
library(shiny)

source("G:/Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/NGS_Rtools_dev/NGS_functions.R")

inputChoices <- loadSeqFolders()
refCohort <- loadRefData()
tumorChoices <- unique(refCohort$metaData$`Tumor type simple`)
tumorChoices <- c("NULL",tumorChoices[order(tumorChoices)])
geneValsAll <- geneValsType <- NULL
sampleChoices <- "First select seq run"

# Define UI for dataset viewer app ----
ui <- navbarPage("PMC NGS R tools", 
                 tabPanel(tags$b("Metadata and reports"),              
                          # row layout with input and output definitions ----
                          fluidRow(
                            column(4,titlePanel(h4("Generate metadata and get files")),selectInput("seqrunMeta", "Choose a seq run:",choices = inputChoices)),
                            column(4,titlePanel(h4("Make reports")),selectInput("seqrunMakeReport", "Choose a seq run:",choices = inputChoices)),
                            column(4,titlePanel(h4("Merge reports")),selectInput("seqrunMergeReport", "Choose a seq run:",choices = inputChoices))
                          ),
                          fluidRow(
                            column(4,selectInput("typeMetadata", "Choose a type:",choices = c("WTS","WES"))),
                            column(4,selectInput("typeMakeReport", "Choose a type:",choices = c("WTS","WES"))),
                            column(4,selectInput("typeMergeReport", "Choose a type:",choices = c("WTS","WES")))
                          ),
                          fluidRow(
                            column(2,actionButton("getMetadata", tags$b("Make metadata"), icon("paper-plane"), style="color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(2,actionButton("refreshFolders", "Refresh")),
                            column(4,actionButton("makeReport", tags$b("Make reports"), icon("paper-plane"), style="color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(4,actionButton("mergeReport", tags$b("Merge reports"), icon("paper-plane"), style="color: #fff; background-color: #fd8723; border-color: #ffffff"))
                          ),
                          fluidRow(
                            br(),
                            br(),
                            mainPanel(width = 12,tableOutput("messages"))
                          )
                 ),
                 tabPanel(tags$b("Check expression"),
                          
                          fluidRow(
                            column(4,titlePanel(h4("Check expression data")))
                          ),
                          fluidRow(
                            column(4,selectInput("seqRunExpression","Choose a seq run",choices = inputChoices)),
                            column(3,selectInput("sampleExpression","Choose a sample",choices = sampleChoices)),
                            column(3,selectInput("tumorType","Choose a tumor type",choices = tumorChoices)),
                            column(2,textInput("geneExpression","Gene",value = NULL))
                          ),
                          fluidRow(
                            column(4,actionButton("checkExpression", tags$b("Plot expression"), style="color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(4,checkboxInput("makePdfExpression","make PDF"))
                          ),
                          fluidRow(
                            tags$head(tags$style('#my_tooltip1 {background-color: rgba(255,255,255,0.8);position: absolute;width: 600px;z-index: 100;}')),
                            tags$script('$(document).ready(function(){
                                        // id of the plot
                                        $("#plotExpression1").mousemove(function(e){ // ID of uiOutput
                                        $("#my_tooltip1").show();
                                        $("#my_tooltip1").css({
                                        top: (e.pageY - 250) + "px",
                                        left: (e.pageX - 250) + "px"
                                        });
                                        });
                                        });
                                        '),

                            br(),
                            br(),
                            mainPanel(
                              width = 12,
                              plotOutput("plotExpression1",hover = hoverOpts(id = "plot_hover1", delay = 0)),
                              uiOutput("my_tooltip1")
                              )
                          ),
                           fluidRow(
                            tags$head(tags$style('#my_tooltip2 {background-color: rgba(255,255,255,0.8);position: absolute;width: 600px;z-index: 100;}')),
                            tags$script('$(document).ready(function(){
                                        // id of the plot
                                        $("#plotExpression2").mousemove(function(e){ // ID of uiOutput
                                        $("#my_tooltip2").show();
                                        $("#my_tooltip2").css({
                                        top: (e.pageY - 700) + "px",
                                        left: (e.pageX - 250) + "px"
                                        });
                                        });
                                        });
                                        '),

                            br(),
                            mainPanel(
                              width = 12,
                              plotOutput("plotExpression2",hover = hoverOpts(id = "plot_hover2", delay = 0)),
                              uiOutput("my_tooltip2")
                            )
                          )
                          
                 ),
                 tabPanel(tags$b("Mutational Signatures"),
                          fluidRow(
                            column(4,titlePanel(h4("Calculate Mutational Signatures")))
                          ),
                          fluidRow(
                            column(5,selectInput("seqRunSignature","Choose a seq run",choices = inputChoices)),
                            column(7,selectInput("sampleSignature","Choose a sample",choices = sampleChoices,width = "150%"))
                          ),
                          fluidRow(
                            column(5,actionButton("getVcfs", tags$b("Get vcfs"))),
                            column(3,actionButton("calculateSignature", tags$b("Calculate Signature"), style="color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(2,checkboxInput("VAF005","VAF >5%")),
                            column(2,checkboxInput("makePdfSignature","make PDF"))
                          ),
                          fluidRow(
                            br(),
                            br(),
                            mainPanel(width = 12,plotOutput("plotSignature")
                            )
                          )
                 ),
                tags$style(HTML(".navbar-default .navbar-brand {color: #ffffff;}
                                 .navbar { background-color: #fd8723;}
                                 .navbar-default .navbar-nav > li > a {color:#ffffff;}
                                 "))
)

#                                  .navbar-default .navbar-nav > .active > a:hover {color: #ffffff;background-color: fca768;}
#                                 .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
#.navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
#.navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
#                                 .navbar-default .navbar-nav > .active > a,
#.navbar-default .navbar-nav > .active > a:focus,


# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  observeEvent(input$getMetadata,ignoreInit = T,{
    type <- input$typeMetadata
    if(type == "WTS"){
      showModal(modalDialog("Making excel sheets and metadata", footer=NULL))
      metaData <- makeFusionExcelFiles(seqRunDir = input$seqrunMeta)
      removeModal()
      if (!is.data.frame(metaData)){
        showNotification("Word document with Data links not found",type="error")
        metaData <- as.data.frame("Word document with Data links not found")
        names(metaData) <- "Error"
      }else{
        showModal(modalDialog("Downloading expression data... this may take some time", footer=NULL))
        checkExpressionData(folder = input$seqrunMeta)
        removeModal()
      }
    } 
    if(type == "WES"){
      showModal(modalDialog("Making metadata sheet and downloading vcfs", footer=NULL))
      metaData <- makeMetaDataWES(seqRunDir = input$seqrunMeta)
      removeModal()
      if (!is.data.frame(metaData)){
        showNotification("Word document with Data links not found",type="error")
        metaData <- as.data.frame("Word document with Data links not found")
        names(metaData) <- "Error"
      }
    }
    output$messages <- renderTable({
      return(metaData)
    })
  })
  
  observeEvent(input$refreshFolders,ignoreInit = T,{
    updateSelectInput(session,"seqrunMeta",choices=loadSeqFolders() )
    updateSelectInput(session,"seqrunMakeReport",choices=loadSeqFolders() )
    updateSelectInput(session,"seqrunMergeReport",choices=loadSeqFolders() )
  })
  observeEvent(input$makeReport,ignoreInit = T,{
    type <- input$typeMakeReport
    showModal(modalDialog("Generating reports", footer=NULL))
    generateReport(folder = input$seqrunMakeReport,type = type)
    removeModal()
    output$messages <- renderTable({
      return("Succesfully made reports")
    })
  })
  
  observeEvent(input$mergeReport,ignoreInit = T,{
     type <- input$typeMergeReport
     showModal(modalDialog("Merging reports", footer=NULL))
     mergeReports(folder = input$seqrunMergeReport,type = type)
     removeModal()
     output$messages <- renderTable({
       return("Succesfully merged reports")
     })
  })
  
  observe({
    seqRun <- input$seqRunExpression
    WTSoverview <- loadRNAseqOverview(folder=seqRun)
    updateSelectInput(session,"sampleExpression",choices=WTSoverview$`Biomaterial ID` )
  })
  
  observe({
    seqRun <- input$seqRunSignature
    sigChoices <- list.files(paste0(baseDirWES,seqRun,"/rawData/"),pattern="vcf.gz")
    if (length(sigChoices) == 0){
      sigChoices <- "Click get vcfs to download vcf files"
    }
    updateSelectInput(session,"sampleSignature",choices=sigChoices)
  })
  
  observeEvent(input$checkExpression,ignoreInit = T,{
    seqRun <- input$seqRunExpression
    sampleName <- input$sampleExpression
    tumorType <- input$tumorType
    gene <- input$geneExpression
    runData <- loadRunExpressionData(folder = seqRun)
    if(input$makePdfExpression){
#      pdf(paste0("G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WTS (RNA-Seq)/",seqRun,"/",sampleName,"_",gene,"_",tumorType,".pdf"),width = 12,height = 7)
      plotExpression(gene = gene,sample = sampleName,tumorType = tumorType,refData = refCohort,runData = runData,pdf=T,folder=seqRun)
#      dev.off()
    }
    #expPlot <- plotExpression(gene = gene,sample = sample,tumorType = tumorType,refData = refCohort,runData = runData)
    output$plotExpression1 <- renderPlot({
      plotExpression(gene = gene,sample = sampleName,tumorType = tumorType,refData = refCohort,runData = runData,showAll=F)
    })
    output$plotExpression2 <- renderPlot({
      plotExpression(gene = gene,sample = sampleName,tumorType = tumorType,refData = refCohort,runData = runData,showAll=T)
    })
    geneValsAll <- cbind(refCohort$metaData,refCohort$counts[gene,])
    colnames(geneValsAll)[ncol(geneValsAll)] <- "gene"
    if (!(sampleName %in% rownames(geneValsAll))){
      geneValsAll <- rbind(geneValsAll,c(NA,NA,NA,runData[gene,sampleName]))
      rownames(geneValsAll)[nrow(geneValsAll)] <- sampleName
    }
    geneValsAll <- geneValsAll[order(geneValsAll$gene),]
    geneValsAll$order <- c(1:nrow(geneValsAll))
    geneValsAll$PMABM <- rownames(geneValsAll)
    geneValsAll <<- geneValsAll
    geneValsType <- geneValsAll[geneValsAll$`Tumor type simple`==tumorType,]
    if (!(sampleName %in% rownames(geneValsType))){
      geneValsType <- rbind(geneValsType,c(NA,NA,NA,runData[gene,sampleName]))
      rownames(geneValsType)[nrow(geneValsType)] <- sampleName
    }
    geneValsType <- geneValsType[order(geneValsType$gene),]
    geneValsType$order <- c(1:nrow(geneValsType))
    geneValsType <<- geneValsType
  })
  
  output$my_tooltip2 <- renderUI({
    hover <- input$plot_hover2 
    y <- nearPoints(geneValsAll, input$plot_hover2, xvar = "gene", yvar = "order",threshold = 10)
    req(nrow(y) != 0)
    #verbatimTextOutput("vals")
    tableOutput("vals2")
  })
  
  output$my_tooltip1 <- renderUI({
    hover <- input$plot_hover1 
    y <- nearPoints(geneValsType, input$plot_hover1, xvar = "gene", yvar = "order",threshold=10)
    req(nrow(y) != 0)
    #verbatimTextOutput("vals")
    tableOutput("vals1")
  })
  
  output$vals2 <- renderTable({
    hover <- input$plot_hover2 
    y <- nearPoints(geneValsAll, input$plot_hover2, xvar = "gene", yvar = "order",threshold=10)
    # y <- nearPoints(data(), input$plot_hover)["wt"]
    req(nrow(y) != 0)
    # y is a data frame and you can freely edit content of the tooltip 
    # with "paste" function
    y <- y[,c("PMABM","Tumor type simple","Resultaat RNA seq (relevante)")]
    colnames(y)[c(2,3)] <- c("Tumortype","Fusie")
    rownames(y) <- NULL
    return((y))
  })
  
  output$vals1 <- renderTable({
    hover <- input$plot_hover1 
    y <- nearPoints(geneValsType, input$plot_hover1, xvar = "gene", yvar = "order",threshold=10)
    # y <- nearPoints(data(), input$plot_hover)["wt"]
    req(nrow(y) != 0)
    # y is a data frame and you can freely edit content of the tooltip 
    # with "paste" function
    y <- y[,c("PMABM","Tumor type simple","Resultaat RNA seq (relevante)")]
    colnames(y)[c(2,3)] <- c("Tumortype","Fusie")
    rownames(y) <- NULL
    return((y))
  })
  
  observeEvent(input$getVcfs,ignoreInit = T,{
    seqRun <- input$seqRunSignature
    showModal(modalDialog("Downloading vcf files", footer=NULL))
    downloadMutect2vcf(seqRun)
    sigChoices <- list.files(paste0(baseDirWES,seqRun,"/rawData/"),pattern="vcf.gz")
    updateSelectInput(session,"sampleSignature",choices=sigChoices)
    removeModal()
  })
  
  observeEvent(input$calculateSignature,ignoreInit = T,{
    VAF005 <- input$VAF005
    seqRun <- input$seqRunSignature
    vcfFile <- input$sampleSignature
    sampleName <- paste(strsplit(vcfFile,"_")[[1]][c(1,2)],collapse = "_")
    showModal(modalDialog("Processing vcf file", footer=NULL))
    mut_mat <- processVcf(folder = seqRun,vcfFile = vcfFile,VAF005 = VAF005)
    if(input$makePdfSignature){
      if(VAF005){
        pdf(paste0(baseDirWES,seqRun,"/",sampleName,"_mutSignatures_VAF005.pdf"),width = 12,height = 7)
      }else{
        pdf(paste0(baseDirWES,seqRun,"/",sampleName,"_mutSignatures.pdf"),width = 12,height = 7)
      }
      getMutationalSignature(mut_mat,sample_names = sampleName,VAF005 = VAF005)
      dev.off()
    }

    output$plotSignature <- renderPlot({
      getMutationalSignature(mut_mat,sample_names = sampleName,VAF005 = VAF005)
    })
    removeModal()
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


