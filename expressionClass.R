library(umap)

refCohort <- loadRefData()

generateUmapData <- function(refCohort,nFeatures=5000,nComp=100){
  baseDir <- paste0(baseDirWTS,"QualityControl/expressionData/")
  
  if(file.exists(paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version)))){
    return(readRDS(paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version))))
  }else{
    ## log tranfsorm ##
    data <- apply(refCohort$counts,2,function(x) (x/sum(x))*1000000)
    dataLog <- log(data+1)
    
    ## select variable features? ##
    varGenes <- apply(dataLog,1,var)
    meanGenes <- apply(dataLog,1,mean)
    
    varFeatures <- names(varGenes)[order(varGenes,decreasing = T)][c(1:2000)]
    
    ## scale data (subtract mean and divide by variance) and store scaling factors ##
    dataScale <- apply(dataLog,2,function(x) (x-meanGenes)/varGenes)
    
    ## perform PCA and store scaling factors ##
    pr <- prcomp(t(dataScale[varFeatures,]))
    ## to transform new sample do : ((log(newSample+1)-meanGenes)/varGenes)[varFeatures,] %*% pr$rotation

    dataUmap <- umap(pr$x[,c(1:100)],alpha=1,gamma=1)
    umapData <- cbind(rownames(dataUmap$layout),dataUmap$layout,refCohort$metaData)
    colnames(umapData) <- c("PMABM","Dim1","Dim2","Tumortype","Fusion","HIX")
    umapAll <- list("umapData"=umapData,"princomp"=pr,"meanGenes"=meanGenes,"varGenes"=varGenes,"varFeatures"=varFeatures,"umapFull"=dataUmap)
    
    saveRDS(umapAll,paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version)))
    
    return(umapAll)  
  }
}

newSampleCoordinates <- function(classData,newSample){
  pr <- ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation
  predict(classData$umapFull,t(pr[c(1:100)]))
}

plotExpressionClass <- function(umapData,classData,tumorTypes,input,scalePoints=1){
  cols <- rainbow(6)
  plot(umapData$Dim1,umapData$Dim2,pch=20,xlab="Dim1",ylab="Dim2",cex=2.5*scalePoints,col='lightgrey')
  points(umapData$Dim1,umapData$Dim2,pch=20,cex=2*scalePoints,col='darkgrey')
  for ( i in 1:length(tumorTypes)){
    points(umapData[umapData$Tumortype == tumorTypes[i],"Dim1"],umapData[umapData$Tumortype == tumorTypes[i],"Dim2"],col=cols[i],pch=20,cex=2*scalePoints)
  }
  if(input$umapSample != "First select seq run"){
    newData <- read.csv(paste0(baseDirWTS,input$umapDir,"/expressionData/",input$umapDir,"_counts.csv"),sep="\t")
    newCoord <- newSampleCoordinates(classData = classData,newSample = newData[,input$umapSample])
    points(newCoord,pch=20,cex=3*scalePoints,col='black')
    legend("bottomright",legend=c(input$umapSample,tumorTypes),col=c("black",cols),pch=20,bty='n')
  }else{
    legend("bottomright",legend=tumorTypes,col=cols,pch=20,bty='n')
  }
}


