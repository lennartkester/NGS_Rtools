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
    
    varFeatures <- names(varGenes)[order(varGenes,decreasing = T)][c(1:5000)]
    
    ## scale data (subtract mean and divide by variance) and store scaling factors ##
    dataScale <- apply(dataLog,2,function(x) (x-meanGenes)/varGenes)
    
    
    ## perform PCA and store scaling factors ##
    pr <- prcomp(t(dataScale[varFeatures,]))
    ## to transform new sample do : ((log(newSample+1)-meanGenes)/varGenes)[varFeatures,] %*% pr$rotation
    
    ## check sdev of principal components ##
    plot(pr$sdev[c(1:100)])
    ## plot first 2 principal components ##
    plot(pr$x[,1],pr$x[,2])
    
    
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






knn(classData$princomp$x[,c(1:100)],pr)
pr <- apply(newData,2,function(newSample) ((log(((newSample/sum(newSample))*1000000)+1)-classData$meanGenes)/classData$varGenes)[classData$varFeatures] %*% classData$princomp$rotation)
# 
# 
# cols <- rep("black",ncol(dataScale))
# cols[refCohort$metaData$`Tumor type simple` == "B-ALL"] <- 'red'
# cols[refCohort$metaData$`Tumor type simple` == "T-ALL"] <- 'blue'
# cols[refCohort$metaData$`Tumor type simple` == "Neuroblastoma"] <- 'green'
# cols[refCohort$metaData$`Tumor type simple` == "Burkitt lymphoma"] <- 'purple'
# cols[refCohort$metaData$`Tumor type simple` == "Pilocytic astrocytoma"] <- 'orange'
# cols[refCohort$metaData$`Tumor type simple` == "AML"] <- 'lightblue'
# cols[refCohort$metaData$`Tumor type simple` == "Hodgkin lymphoma"] <- 'pink'
# 
# 
# pdf("20200615_classPlotExample.pdf")
# plot(dataUmap$layout,col=cols,pch=20)
# legend("bottomleft",pch=20,bty='n',legend=c("B-ALL","T-ALL","Neuroblastoma","Burkitt lymphoma","Pilocytic astrocytoma","AML","Hodgkin lymphoma"),col=c('red','blue','green','purple','orange','lightblue','pink'))
# dev.off()

