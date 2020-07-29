
library(shiny)

scriptDir <- "G:/Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/NGS_Rtools_dev/"

source(paste0(scriptDir,"NGS_functions_general.R"))
source(paste0(scriptDir,"NGS_functions_metadata.R"))
source(paste0(scriptDir,"NGS_functions_reports.R"))
source(paste0(scriptDir,"NGS_functions_expression.R"))
source(paste0(scriptDir,"NGS_functions_mutSignatures.R"))
source(paste0(scriptDir,"NGS_functions_RNAclassifier.R"))

inputChoices <- loadSeqFolders()
refCohort <- loadRefData()
classData <- generateUmapData(refCohort = refCohort)
umapData <- classData$umapData
tumorChoices <- unique(refCohort$metaData$Disease_sub_class)
tumorChoices <- c("NULL",tumorChoices[order(tumorChoices)])
subTypeChoices <- unique(refCohort$metaData$Disease_sub_specification1)
subTypeChoices <- c("NULL",subTypeChoices[order(subTypeChoices)])

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
                            column(4,selectInput("typeMakeReport", "Choose a type:",choices = c("WTS","WES")))
                            #column(4,selectInput("typeMergeReport", "Choose a type:",choices = c("WTS","WES")))
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
                 tabPanel(tags$b("Compare expression"),
                          fluidRow(
                            column(4,titlePanel(h4("Highlight tumor types")))
                          ),
                          fluidRow(
                            column(2,selectInput("tumorType1", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType1", "Choose a Subtype:",choices = subTypeChoices)),
                            column(2,selectInput("tumorType2", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType2", "Choose a Subtype:",choices = subTypeChoices)),
                            column(2,selectInput("tumorType3", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType3", "Choose a Subtype:",choices = subTypeChoices))
                          ),
                          fluidRow(
                            column(2,selectInput("tumorType4", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType4", "Choose a Subtype:",choices = subTypeChoices)),
                            column(2,selectInput("tumorType5", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType5", "Choose a Subtype:",choices = subTypeChoices)),
                            column(2,selectInput("tumorType6", "Choose a Tumortype:",choices = tumorChoices)),
                            column(2,selectInput("tumorSubType6", "Choose a Subtype:",choices = subTypeChoices))
                          ),
                          fluidRow(
                            column(2,actionButton("resetClassPlot",tags$b("Reset"), style="margin-top: 25px;color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(2,actionButton("highlightTypes",tags$b("Highlight tumortypes"), style="margin-top: 25px;color: #fff; background-color: #fd8723; border-color: #ffffff"))
                          ),
                          fluidRow(
                            column(4,titlePanel(h4("Highlight gene")))
                          ),
                          fluidRow(
                            column(3,textInput("geneExpressionClass","Gene",value = NULL)),
                            column(2,actionButton("highlightGeneClass",tags$b("Highlight gene"), style="margin-top: 25px;color: #fff; background-color: #fd8723; border-color: #ffffff"))
                          ),
                          fluidRow(
                            column(4,titlePanel(h4("Add new sample")))
                          ),
                          fluidRow(
                            column(3,selectInput("umapDir","Choose a seq run",choices = inputChoices)),
                            column(3,selectInput("umapSample","Choose a sample",choices = sampleChoices)),
                            column(3,actionButton("addSampleExpClass",tags$b("Classify sample"), style="margin-top: 25px;color: #fff; background-color: #fd8723; border-color: #ffffff")),
                            column(3,actionButton("printExpClass",tags$b("Print PDF"), style="margin-top: 25px;color: #fff; background-color: #fd8723; border-color: #ffffff"))
                           ),
                          fluidRow(
                            mainPanel(width = 6,tableOutput("classTable2")),
                            mainPanel(width = 6,tableOutput("classTable"))
                          ),
                          fluidRow(
                            style = "height:600px",
                            tags$head(tags$style('#my_tooltip3 {background-color: rgba(255,255,255,0.8);position: absolute;width: 600px;z-index: 100;}')),
                            tags$script('$(document).ready(function(){
                                        // id of the plot
                                        $("#compareExpression").mousemove(function(e){ // ID of uiOutput
                                          $("#my_tooltip3").show();
                                          $("#my_tooltip3").css({
                                          top: (e.pageY - 600) + "px",
                                          left: (e.pageX - 50) + "px"
                                          })

                                        });

                                       });
                                      
                                        '),
                            
                            br(),
                            mainPanel(
                              width = 12,
                              plotOutput("compareExpression",hover = hoverOpts(id = "plot_hover3", delay = 0),inline=T),
                              uiOutput("my_tooltip3")
                            ),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            br()
                            
                            
                          )
                 ),
                tags$style(HTML(".navbar-default .navbar-brand {color: #ffffff;}
                                 .navbar { background-color: #fd8723;}
                                 .navbar-default .navbar-nav > li > a {color:#ffffff;}
                                 ")),
                tags$head(
                  tags$style(
                    HTML(".shiny-notification {
                         position:fixed;
                         width:300px;  
                         top: calc(22px);
                         left: calc(250px);
                         opacity: 1;
                         }
                         "
                    )
                  )
                )
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
    out <- generateReport(folder = input$seqrunMakeReport,type = type)
    removeModal()
    if(out != "Succesfully generated reports"){
      showNotification(paste(out),duration = 10,closeButton = T,type = "error")
    }else{
      showNotification(paste(out),duration = 5,closeButton = T,type = "default")
      
    }
    #output$messages <- renderTable({
    #  return(out)
    #})
  })
  
  observeEvent(input$mergeReport,ignoreInit = T,{
     showModal(modalDialog("Merging reports", footer=NULL))
     out <- tryCatch(mergeReports(folder = input$seqrunMergeReport),error=function(e) return(paste("Could not combine reports, destination file open?")))
     removeModal()
     
     if(out != "Succesfully merged reports"){
       showNotification(paste(out),duration = 10,closeButton = T,type = "error")
     }else{
       showNotification(paste(out),duration = 5,closeButton = T,type = "default")
       
     }
  })
  
  observe({
    seqRun <- input$seqRunExpression
    WTSoverview <- loadRNAseqOverview(folder=seqRun)
    updateSelectInput(session,"sampleExpression",choices=WTSoverview$`Biomaterial ID` )
  })
  
  observe({
    seqRun <- input$seqRunSignature
    sigChoices <- list.files(paste0(baseDirWES,seqRun,"/rawData/"),pattern="vcf.gz")
    sigChoices <- sigChoices[sapply(sigChoices,function(x) length(strsplit(x,"_")[[1]])) == 4]
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
  
  output$my_tooltip1 <- renderUI({
    hover <- input$plot_hover1 
    y <- nearPoints(geneValsType, input$plot_hover1, xvar = "gene", yvar = "order",threshold=10)
    req(nrow(y) != 0)
    #verbatimTextOutput("vals")
    tableOutput("vals1")
  })
  
  output$my_tooltip2 <- renderUI({
    hover <- input$plot_hover2 
    y <- nearPoints(geneValsAll, input$plot_hover2, xvar = "gene", yvar = "order",threshold = 10)
    req(nrow(y) != 0)
    #verbatimTextOutput("vals")
    tableOutput("vals2")
  })
  
  output$vals1 <- renderTable({
    hover <- input$plot_hover1 
    y <- nearPoints(geneValsType, input$plot_hover1, xvar = "gene", yvar = "order",threshold=10)
    # y <- nearPoints(data(), input$plot_hover)["wt"]
    req(nrow(y) != 0)
    # y is a data frame and you can freely edit content of the tooltip 
    # with "paste" function
    y <- y[,c("PMABM","Disease_sub_specification1","Resultaat RNA seq (relevante)")]
    colnames(y)[c(2,3)] <- c("Tumortype","Fusie")
    rownames(y) <- NULL
    return((y))
  })
  
  output$vals2 <- renderTable({
    hover <- input$plot_hover2 
    y <- nearPoints(geneValsAll, input$plot_hover2, xvar = "gene", yvar = "order",threshold=10)
    # y <- nearPoints(data(), input$plot_hover)["wt"]
    req(nrow(y) != 0)
    # y is a data frame and you can freely edit content of the tooltip 
    # with "paste" function
    y <- y[,c("PMABM","Disease_sub_specification1","Resultaat RNA seq (relevante)")]
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
  
  classPlotInput <- reactiveValues(tumorTypes=NULL,classSample=NULL,geneExpressionClass=NULL,neighbours=NULL)
  
  observeEvent(input$resetClassPlot,ignoreInit = T,{
    classPlotInput$tumorTypes <- NULL
    classPlotInput$classSample <- NULL
    classPlotInput$geneExpressionClass <- NULL
    classPlotInput$neighbours <- NULL
    
  })
  
  observeEvent(input$highlightTypes, ignoreInit = T, {
    tumorTypes <- matrix(ncol=2,nrow=6,data=NA)
    tumorTypes[1,] <- c(input$tumorType1,input$tumorSubType1)
    tumorTypes[2,] <- c(input$tumorType2,input$tumorSubType2)
    tumorTypes[3,] <- c(input$tumorType3,input$tumorSubType3)
    tumorTypes[4,] <- c(input$tumorType4,input$tumorSubType4)
    tumorTypes[5,] <- c(input$tumorType5,input$tumorSubType5)
    tumorTypes[6,] <- c(input$tumorType6,input$tumorSubType6)
    classPlotInput$tumorTypes <- tumorTypes
    classPlotInput$geneExpressionClass <- NULL
  })
  
  observeEvent(input$highlightGeneClass, ignoreInit = T, {
    classPlotInput$geneExpressionClass <- input$geneExpressionClass
    classPlotInput$tumorTypes <- NULL
    classPlotInput$classSample <- NULL
    classPlotInput$neighbours <- NULL
  })
  
  observeEvent(input$addSampleExpClass, ignoreInit = T, {
    neighbours <- NULL
    if(input$umapSample != "First select seq run"){
      classResults <- predictClass(input = input,classData = classData)
      res <- classResults$results
      res <- as.matrix(res[,order(res,decreasing = T)[c(1:3)]])
      colnames(res) <- input$umapSample
      html_caption_str <- as.character(shiny::tags$b(style = "color: #fd8723", "Predicted tumortype sub specification"))
      output$classTable <- renderTable(rownames = T,{
        t(res)/100
      }, caption = html_caption_str, caption.placement = "top")
      res2 <- classResults$results2
      res2 <- as.matrix(res2[,order(res2,decreasing = T)[c(1:3)]])
      colnames(res2) <- input$umapSample
      html_caption_str <- as.character(shiny::tags$b(style = "color: #fd8723", "Predicted tumortype"))
      output$classTable2 <- renderTable(rownames = T,{
        t(res2)/100
      }, caption = html_caption_str, caption.placement = "top")
      classPlotInput$neighbours <- classResults$neighbours
      classPlotInput$classSample <- list(umapDir=input$umapDir,umapSample=input$umapSample)
      classPlotInput$geneExpressionClass <- NULL
    }    
  })

  output$compareExpression <- renderPlot(height=900,width=1200,{
    par(mar=c(10,4,3,3))
    plotExpressionClass(umapData = umapData,
                        classData = classData,
                        tumorTypes = classPlotInput$tumorTypes,
                        umapDir=classPlotInput$classSample$umapDir,
                        umapSample=classPlotInput$classSample$umapSample,
                        neighbours=classPlotInput$neighbours,
                        geneExpressionClass=classPlotInput$geneExpressionClass)
  })
  
  observeEvent(input$printExpClass,ignoreInit=T, {
    tumorTypes <- c(input$tumorType1,input$tumorType2,input$tumorType3,input$tumorType4,input$tumorType5,input$tumorType6)
    if(input$umapSample != "First select seq run"){
      pdf(paste0(baseDirWTS,input$umapDir,"/",input$umapSample,"_expressionClass.pdf"),width = 12,height = 7)
      plotExpressionClass(umapData = umapData,classData = classData,tumorTypes = tumorTypes,input=input,scalePoints=0.6)
      dev.off()
    }else{
      showNotification("First select seq run and sample",duration = 5,type="error")
    }
  })
  
  output$vals3 <- renderTable({
    hover <- input$plot_hover3 
    y <- nearPoints(umapData, input$plot_hover3, xvar = "Dim1", yvar = "Dim2",threshold=5)
    # y <- nearPoints(data(), input$plot_hover)["wt"]
    req(nrow(y) != 0)
    # y is a data frame and you can freely edit content of the tooltip 
    # with "paste" function
    y <- y[,c("PMABM","Disease_sub_specification1","Fusion","HIX")]
    if (nrow(y) > 5){
      y <- y[c(1:5),]
    }
    rownames(y) <- NULL
    return((y))
  })
  
  output$my_tooltip3 <- renderUI({
    hover <- input$plot_hover3 
    y <- nearPoints(umapData, input$plot_hover3, xvar = "Dim1", yvar = "Dim2",threshold = 5)
    req(nrow(y) != 0)
    tableOutput("vals3")
  })
  
 observe({
   seqRun <- input$umapDir
   WTSoverview <- loadRNAseqOverview(folder=seqRun)
   updateSelectInput(session,"umapSample",choices=WTSoverview$`Biomaterial ID` )
 })
 
 observe({
   tumorType <- input$tumorType1
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType1",choices=subChoices )
 })
 observe({
   tumorType <- input$tumorType2
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType2",choices=subChoices )
 })
 observe({
   tumorType <- input$tumorType3
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType3",choices=subChoices )
 })
 observe({
   tumorType <- input$tumorType4
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType4",choices=subChoices )
 })
 observe({
   tumorType <- input$tumorType5
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType5",choices=subChoices )
 })
 observe({
   tumorType <- input$tumorType6
   subChoices <- unique(refCohort$metaData$Disease_sub_specification1[refCohort$metaData$Disease_sub_class == tumorType])
   subChoices <- subChoices[order(subChoices)]
   if (length(subChoices) > 1){
     subChoices <- c("All subtypes",subChoices)
   }
   updateSelectInput(session,"tumorSubType6",choices=subChoices )
 })
   
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


