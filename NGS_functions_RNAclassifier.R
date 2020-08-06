library(umap)

refCohort <- loadRefData()

generateUmapData <- function(refCohort,nFeatures=5000,nComp=100,domain=T){
  baseDir <- paste0(baseDirWTS,"QualityControl/expressionData/")
  if(domain){
    domains <- c("All","Neuro","Solid","Hemato")
    umapAll <- list()
    for ( d in domains){
      if(file.exists(paste0(baseDir,sub(".csv",paste0("_umapData_",d,".rds"),refCohort$version)))){
        umapAll[[d]] <- readRDS(paste0(baseDir,sub(".csv",paste0("_umapData_",d,".rds"),refCohort$version)))
      }else{
        ## log tranfsorm ##
        if (d == "All"){
          data <- refCohort$counts
          metaData <- refCohort$metaData
        }else{
          data <- refCohort$counts[,refCohort$metaData$Domain == d]
          metaData <- refCohort$metaData[refCohort$metaData$Domain == d,]
        }
        
        data <- apply(data,2,function(x) (x/sum(x))*1000000)
        dataLog <- log(data+1)
        
        ## select variable features ##
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
        umapData <- cbind(rownames(dataUmap$layout),dataUmap$layout,metaData)
        colnames(umapData)[c(1:6)] <- c("PMABM","Dim1","Dim2","Tumortype","Fusion","HIX")
        umapDomain <- list("umapData"=umapData,"princomp"=pr,"meanGenes"=meanGenes,"varGenes"=varGenes,"varFeatures"=varFeatures,"umapFull"=dataUmap)
        
        umapData$Tumortype <- as.factor(umapData$Tumortype)
        umapData$Disease_sub_specification1 <- as.factor(umapData$Disease_sub_specification1)
        
        kMat <- as.data.frame(matrix(ncol=2,nrow=100))
        samplesTrainAll <- list()
        k <- list()
        for ( j in c(1:100)){
          #  set.seed(j)
          #  samplesTrain <- unique(sample(rownames(prData),replace = T))
          samplesTrain <- rownames(prData)
          set.seed(j)
          compsTrain <- unique(sample(c(1:nComp),replace = T))
          samplesTrainAll[[j]] <- samplesTrain
          prdataTrain <- as.data.frame(prData[samplesTrain,compsTrain])
          prdataTrain$class <- umapData[samplesTrain,"Disease_sub_specification1"]
          #prdataTest <- as.data.frame(pr$x[samplesTest,c(1:100)])
          #names(prCurSample) <- colnames(prdataTrain)[compsTrain]
          #k[[j]] <- kknn(class~., prdata,as.data.frame(prCurSample),distance = 2,kernel = "optimal",scale=F,k=10)
          kTrain <- train.kknn(class~., prdataTrain,distance = 2,kernel = "optimal",scale=F,kmax = 20)
          k[[j]] <- kTrain
          kMat[j,1] <- kTrain$best.parameters$k
          kMat[j,2] <- kTrain$best.parameters$kernel
          
        }
        trainFreqs <- table(unlist(samplesTrainAll))
        umapDomain$kMat <- kMat
        umapDomain$trainFreqs <- trainFreqs
        
        saveRDS(umapDomain,paste0(baseDir,sub(".csv",paste0("_umapData_",d,".rds"),refCohort$version)))
        umapAll[[d]] <- umapDomain
      }
    }
    return(umapAll)  
  }else{
    if(file.exists(paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version)))){
      return(readRDS(paste0(baseDir,sub(".csv","_umapData.rds",refCohort$version))))
    }else{
      ## log tranfsorm ##
      data <- apply(refCohort$counts,2,function(x) (x/sum(x))*1000000)
      dataLog <- log(data+1)
      
      ## select variable features ##
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
      colnames(umapData)[c(1:6)] <- c("PMABM","Dim1","Dim2","Tumortype","Fusion","HIX")
      umapAll <- list("umapData"=umapData,"princomp"=pr,"meanGenes"=meanGenes,"varGenes"=varGenes,"varFeatures"=varFeatures,"umapFull"=dataUmap)
      
      umapData$Tumortype <- as.factor(umapData$Tumortype)
      umapData$Disease_sub_specification1 <- as.factor(umapData$Disease_sub_specification1)
      
      kMat <- as.data.frame(matrix(ncol=2,nrow=100))
      samplesTrainAll <- list()
      k <- list()
      for ( j in c(1:100)){
        #  set.seed(j)
        #  samplesTrain <- unique(sample(rownames(prData),replace = T))
        samplesTrain <- rownames(prData)
        set.seed(j)
        compsTrain <- unique(sample(c(1:nComp),replace = T))
        samplesTrainAll[[j]] <- samplesTrain
        prdataTrain <- as.data.frame(prData[samplesTrain,compsTrain])
        prdataTrain$class <- umapData[samplesTrain,"Disease_sub_specification1"]
        #prdataTest <- as.data.frame(pr$x[samplesTest,c(1:100)])
        #names(prCurSample) <- colnames(prdataTrain)[compsTrain]
        #k[[j]] <- kknn(class~., prdata,as.data.frame(prCurSample),distance = 2,kernel = "optimal",scale=F,k=10)
        kTrain <- train.kknn(class~., prdataTrain,distance = 2,kernel = "optimal",scale=F,kmax = 20)
        k[[j]] <- kTrain
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
}

newSampleCoordinates <- function(classData,newSample){
  pr <- ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation
  predict(classData$umapFull,t(pr[c(1:100)]))
}

newSamplePrinComp <- function(classData,newSample){
  ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation
}

plotExpressionClass <- function(classData,tumorTypes,umapDir=NULL,umapSample=NULL,scalePoints=1,neighbours=NULL,geneExpressionClass=NULL){
  umapData <- classData$umapData
  cols <- rainbow(6)
  ymin <- min(umapData$Dim2)-((max(umapData$Dim2)-min(umapData$Dim2))*0.05)
  ymax <- max(umapData$Dim2)+((max(umapData$Dim2)-min(umapData$Dim2))*0.05)
  xmin <- min(umapData$Dim1)-((max(umapData$Dim1)-min(umapData$Dim1))*0.05)
  xmax <- max(umapData$Dim1)+((max(umapData$Dim1)-min(umapData$Dim1))*0.05)
  plot(umapData$Dim1,umapData$Dim2,pch=20,xlab="Dim1",ylab="Dim2",ylim=c(ymin,ymax),xlim=c(xmin,xmax),cex=2.5*scalePoints,col='lightgrey')
  points(umapData$Dim1,umapData$Dim2,pch=20,cex=2*scalePoints,col='darkgrey')
  legendMatrix <- matrix(ncol=2,nrow=6,data=NA)
  if(!is.null(tumorTypes)){
    for ( i in 1:nrow(tumorTypes)){
      if(tumorTypes[i,1] == "Select type"){
        next
      }
      if(tumorTypes[i,2] == "All subtypes"){
        points(umapData[umapData$Disease_sub_class == tumorTypes[i,1],"Dim1"],umapData[umapData$Disease_sub_class == tumorTypes[i,1],"Dim2"],col=cols[i],pch=20,cex=2*scalePoints)
        legendMatrix[i,1] <- tumorTypes[i,1]
        legendMatrix[i,2] <- cols[i]
      }else{
        points(umapData[umapData$Disease_sub_specification1 == tumorTypes[i,2],"Dim1"],umapData[umapData$Disease_sub_specification1 == tumorTypes[i,2],"Dim2"],col=cols[i],pch=20,cex=2*scalePoints)
        legendMatrix[i,1] <- tumorTypes[i,2]
        legendMatrix[i,2] <- cols[i]
      }
    }  
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
    legend("bottomright",legend=c(umapSample,legendMatrix[,1]),col=c("black",legendMatrix[,2]),pch=20,bty='n')
  }else{
    if(!is.null(tumorTypes)){
      legend("bottomright",legend=legendMatrix[,1],col=legendMatrix[,2],pch=20,bty='n')
    }
  }
  if(!is.null(neighbours)){
    neighbourPoints <- umapData[names(neighbours)[neighbours > 0],]
    neighbourPoints$neighbour <- neighbours[neighbours > 0]
    neighbourPoints <- neighbourPoints[order(neighbourPoints$neighbour),]
    cols <- rev(terrain.colors(101))
    neighbourPoints$col <- cols[round(neighbourPoints$neighbour*100,0)+1]
    points(neighbourPoints$Dim1,neighbourPoints$Dim2,col=neighbourPoints$col,pch=20,cex=2)
    legendRange <- seq((mean(c(xmin,xmax))-((xmax-xmin)*0.1)),(mean(c(xmin,xmax))+((xmax-xmin)*0.1)),length.out = 101)
    #legendRange <- seq((min(umapData$Dim1)+max(umapData$Dim1))-((max(umapData$Dim1)-min(umapData$Dim1))*0.1),(min(umapData$Dim1)+max(umapData$Dim1))+((max(umapData$Dim1)-min(umapData$Dim1))*0.1),length.out = 101)
    legendy <- ymin
    for(i in 1:100){
      lines(x = legendRange[c(i,(i+1))],rep(legendy,2),lwd=10,col=rev(terrain.colors(100))[i])
    }
    text(x=legendRange[1],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "0")
    text(x=legendRange[101],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "1")
    text(x=legendRange[50],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "Similarity")
  }
  if(!is.null(geneExpressionClass)){
    if(!(geneExpressionClass %in% rownames(refCohort$counts))){
      return()
    }
    geneData <- refCohort$counts[geneExpressionClass,]
    #geneData <- (geneData-mean(geneData))/sd(geneData)
    colorSeq <- c(min(geneData)-0.1,seq(quantile(geneData,c(0.10,0.99))[1],quantile(geneData,c(0.10,0.99))[2],length.out = 98),max(geneData))
    cols <- rev(magma(100))[as.numeric(cut(geneData,colorSeq))]
    names(cols) <- colnames(refCohort$counts)
    points(umapData$Dim1,umapData$Dim2,pch=20,cex=2*scalePoints,col=cols[umapData$PMABM])
    legendRange <- seq((mean(c(xmin,xmax))-((xmax-xmin)*0.1)),(mean(c(xmin,xmax))+((xmax-xmin)*0.1)),length.out = 101)
    #legendRange <- seq((min(umapData$Dim1)+max(umapData$Dim1))-((max(umapData$Dim1)-min(umapData$Dim1))*0.1),(min(umapData$Dim1)+max(umapData$Dim1))+((max(umapData$Dim1)-min(umapData$Dim1))*0.1),length.out = 101)
    legendy <- ymin
    for(i in 1:100){
      lines(x = legendRange[c(i,(i+1))],rep(legendy,2),lwd=10,col=rev(magma(100))[i])
    }
    text(x=legendRange[1],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = paste(round(quantile(geneData,c(0.10,0.99))[1],2)))
    text(x=legendRange[101],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = paste(round(quantile(geneData,c(0.10,0.99))[2],2)))
    text(x=legendRange[50],y=legendy+((max(umapData$Dim2)-min(umapData$Dim2))*0.025),labels = "Counts per million")
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
    #set.seed(j)
    #samplesTrain <- unique(sample(rownames(prData),replace = T))
    samplesTrain <- rownames(prData)
    set.seed(j)
    compsTrain <- unique(sample(c(1:100),replace = T))
    
    #samplesTrainAll[[j]] <- samplesTrain
    prdataTrain <- as.data.frame(prData[samplesTrain,compsTrain])
    prdataTrain$class <- as.factor(classData$umapData[samplesTrain,"Disease_sub_specification1"])
    #prdataTest <- as.data.frame(pr$x[samplesTest,c(1:100)])
    colnames(prCurSample) <- colnames(prData)
    prCurSample2 <- prCurSample[,compsTrain]
    #k[[j]] <- kknn(class~., prdata,as.data.frame(prCurSample),distance = 2,kernel = "optimal",scale=F,k=10)
    #kTrain <- train.kknn(class~., prdataTrain,distance = 2,kernel = "optimal",scale=F,kmax = 20)
    #print(paste(kTrain$best.parameters$k,kTrain$best.parameters$kernel))
    k[[j]] <- kknn(class~., prdataTrain,(as.data.frame(t(prCurSample2))),distance = 2,kernel = as.character(classData$kMat$V2[j]),scale=F,k=classData$kMat$V1[j])
    rownames(k[[j]]$prob) <- input$umapSample
    nb <- as.vector(k[[j]]$W)
    names(nb) <- rownames(prdataTrain)[k[[j]]$C]
    neighbourMat[j,names(nb)] <- nb
    #neighbours[[j]] <- nb
  }
  res <- Reduce('+',lapply(k,function(x) x$prob))
  
  classMatrix <- classData$umapData[,c(7:10)]
  classMatrix <- classMatrix[!duplicated(apply(classMatrix,1,function(x) paste(x,collapse="_"))),]
  
  res2 <- matrix(nrow=1,ncol=length(unique(classData$umapData$Disease_sub_class)),data=0)
  colnames(res2) <- unique(classData$umapData$Disease_sub_class)
  for ( j in 1:ncol(res)){
    subClass <- classMatrix[colnames(res)[j] == classMatrix$Disease_sub_specification1,3]
    res2[,subClass] <- res2[,subClass] + res[,j]   
  }

  #neighbourFreqs <- table(unlist(neighbours))/classData$trainFreqs[names(table(unlist(neighbours)))]
  neighbourScores <- apply(neighbourMat,2,sum,na.rm=T)
  neighbourScores <- neighbourScores/max(neighbourScores)
  return(list("results"=res,"results2"=res2,"neighbours"=neighbourScores))
}



