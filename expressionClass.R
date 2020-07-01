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
    
    varFeatures <- names(varGenes)[order(varGenes,decreasing = T)][c(1:nFeatures)]
    
    ## scale data (subtract mean and divide by variance) and store scaling factors ##
    dataScale <- apply(dataLog,2,function(x) (x-meanGenes)/varGenes)
    
    ## perform PCA and store scaling factors ##
    pr <- prcomp(t(dataScale[varFeatures,]))
    prData <- pr$x
    ## to transform new sample do : ((log(newSample+1)-meanGenes)/varGenes)[varFeatures,] %*% pr$rotation

    dataUmap <- umap(pr$x[,c(1:nComp)],alpha=1,gamma=1)
    umapData <- cbind(rownames(dataUmap$layout),dataUmap$layout,refCohort$metaData)
    colnames(umapData) <- c("PMABM","Dim1","Dim2","Tumortype","Fusion","HIX")
    umapAll <- list("umapData"=umapData,"princomp"=pr,"meanGenes"=meanGenes,"varGenes"=varGenes,"varFeatures"=varFeatures,"umapFull"=dataUmap)
    
    umapData$Tumortype <- as.factor(umapData$Tumortype)
    
    kMat <- as.data.frame(matrix(ncol=2,nrow=100))
    samplesTrainAll <- list()
    for ( j in c(1:100)){
      set.seed(j)
      samplesTrain <- unique(sample(rownames(prData),replace = T))
      samplesTrainAll[[j]] <- samplesTrain
      prdataTrain <- as.data.frame(prData[samplesTrain,c(1:100)])
      prdataTrain$class <- umapData[samplesTrain,"Tumortype"]
      #prdataTest <- as.data.frame(pr$x[samplesTest,c(1:100)])
      names(prCurSample) <- colnames(prdataTrain)[c(1:100)]
      #k[[j]] <- kknn(class~., prdata,as.data.frame(prCurSample),distance = 2,kernel = "optimal",scale=F,k=10)
      kTrain <- train.kknn(class~., prdataTrain,distance = 2,kernel = "optimal",scale=F,kmax = 20)
      kMat[j,1] <- kTrain$best.parameters$k
      kMat[j,2] <- kTrain$best.parameters$kernel
      
    }
    trainFreqs <- table(unlist(samplesTrainAll))
    umapAll$kMat <- kMat
    umapAll$trainFreqs <- trainFreqs
    
    saveRDS(umapAll,paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version)))
    
    return(umapAll)  
  }
}

newSampleCoordinates <- function(classData,newSample){
  pr <- ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation
  predict(classData$umapFull,t(pr[c(1:100)]))
}

newSamplePrinComp <- function(classData,newSample){
  ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation
}

plotExpressionClass <- function(umapData,classData,tumorTypes,umapDir=NULL,umapSample=NULL,scalePoints=1,neighbours=NULL){
  cols <- rainbow(6)
  plot(umapData$Dim1,umapData$Dim2,pch=20,xlab="Dim1",ylab="Dim2",cex=2.5*scalePoints,col='lightgrey')
  points(umapData$Dim1,umapData$Dim2,pch=20,cex=2*scalePoints,col='darkgrey')
  for ( i in 1:length(tumorTypes)){
    points(umapData[umapData$Tumortype == tumorTypes[i],"Dim1"],umapData[umapData$Tumortype == tumorTypes[i],"Dim2"],col=cols[i],pch=20,cex=2*scalePoints)
  }
  if(!is.null(umapSample)){
    if(umapSample %in% umapData$PMABM){
      points(umapData[umapData$PMABM == umapSample,"Dim1"],umapData[umapData$PMABM == umapSample,"Dim2"],pch=20,cex=3*scalePoints,col='black')
    }else{
      checkExpressionData(umapDir)
      newData <- read.csv(paste0(baseDirWTS,umapDir,"/expressionData/",umapDir,"_counts.csv"),sep="\t")
      newCoord <- newSampleCoordinates(classData = classData,newSample = newData[,umapSample])
      points(newCoord,pch=20,cex=3*scalePoints,col='black')
    }
    legend("bottomright",legend=c(umapSample,tumorTypes),col=c("black",cols),pch=20,bty='n')
  }else{
    if(!is.null(tumorTypes)){
      legend("bottomright",legend=tumorTypes,col=cols,pch=20,bty='n')
    }
  }
  if(!is.null(neighbours)){
    neighbourPoints <- umapData[names(neighbours)[neighbours > 0],]
    neighbourPoints$neighbour <- neighbours[neighbours > 0]
    neighbourPoints <- neighbourPoints[order(neighbourPoints$neighbour),]
    cols <- rev(terrain.colors(101))
    neighbourPoints$col <- cols[round(neighbourPoints$neighbour*100,0)+1]
    points(neighbourPoints$Dim1,neighbourPoints$Dim2,col=neighbourPoints$col,pch=20,cex=2)
    legendRange <- seq((min(umapData$Dim1)+max(umapData$Dim1))-((max(umapData$Dim1)-min(umapData$Dim1))*0.1),(min(umapData$Dim1)+max(umapData$Dim1))+((max(umapData$Dim1)-min(umapData$Dim1))*0.1),length.out = 101)
    legendy <- min(umapData$Dim2)
    for(i in 1:100){
      lines(x = legendRange[c(i,(i+1))],rep(legendy,2),lwd=10,col=rev(terrain.colors(100))[i])
    }
    text(x=legendRange[1],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "0")
    text(x=legendRange[101],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "1")
    text(x=legendRange[50],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "Similarity")
  }
}

predictClass <- function(input,classData){
  checkExpressionData(input$umapDir)
  newData <- read.csv(paste0(baseDirWTS,input$umapDir,"/expressionData/",input$umapDir,"_counts.csv"),sep="\t")
  prCurSample <- newSamplePrinComp(classData = classData,newSample = newData[,input$umapSample])
  prData <- classData$princomp$x
  if(input$umapSample %in% rownames(prData)){
    prData <- prData[rownames(prData) != input$umapSample,]
  }
  classData$umapData$Tumortype <- as.factor(classData$umapData$Tumortype)
  k <- list()
  #samplesTrainAll <- list()
  neighbourMat <- matrix(ncol=nrow(prData),nrow=100)
  colnames(neighbourMat) <- rownames(prData)
  for ( j in c(1:100)){
    set.seed(j)
    samplesTrain <- unique(sample(rownames(prData),replace = T))
    #samplesTrainAll[[j]] <- samplesTrain
    prdataTrain <- as.data.frame(prData[samplesTrain,c(1:100)])
    prdataTrain$class <- classData$umapData[samplesTrain,"Tumortype"]
    #prdataTest <- as.data.frame(pr$x[samplesTest,c(1:100)])
    names(prCurSample) <- colnames(prdataTrain)[c(1:100)]
    #k[[j]] <- kknn(class~., prdata,as.data.frame(prCurSample),distance = 2,kernel = "optimal",scale=F,k=10)
    #kTrain <- train.kknn(class~., prdataTrain,distance = 2,kernel = "optimal",scale=F,kmax = 20)
    #print(paste(kTrain$best.parameters$k,kTrain$best.parameters$kernel))
    k[[j]] <- kknn(class~., prdataTrain,as.data.frame(prCurSample),distance = 2,kernel = as.character(classData$kMat$V2[j]),scale=F,k=classData$kMat$V1[j])
    rownames(k[[j]]$prob) <- input$umapSample
    nb <- as.vector(k[[j]]$W)
    names(nb) <- rownames(prdataTrain)[k[[j]]$C]
    neighbourMat[j,names(nb)] <- nb
    #neighbours[[j]] <- nb
  }
  res <- Reduce('+',lapply(k,function(x) x$prob))
  neighbourFreqs <- table(unlist(neighbours))/classData$trainFreqs[names(table(unlist(neighbours)))]
  neighbourScores <- apply(neighbourMat,2,sum,na.rm=T)
  neighbourScores <- neighbourScores/max(neighbourScores)
  return(list("results"=res,"neighbours"=neighbourScores))
}



