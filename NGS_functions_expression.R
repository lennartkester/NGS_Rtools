

plotExpression <- function(gene=NULL,sample=NULL,tumorType=NULL,folder=NULL,refData=refCohort,runData=runCohort,pdf=F,showAll=T){
  if(!gene %in% rownames(refData$counts)){
    stop(paste(gene,"does not exists in refData, are you sure spelling is correct?"))
  }
  if(!is.null(tumorType)){
    if(!(tumorType %in% refData$metaData$`Tumor type simple`)){
      tumorType2 <- agrep(tumorType,unique(refData$metaData$`Tumor type simple`),value=TRUE,ignore.case = TRUE)[1]
      if(is.na(tumorType2)){
        stop(cat(paste("tumorType could not be matched\ntumorType can be one of the following:\n",paste(unique(refData$metaData$`Tumor type simple`[order(refData$metaData$`Tumor type simple`)]),collapse="\t\t"),sep="")))
      }
      message(paste0("Could not find ",tumorType,", assuming you mean ",tumorType2))
      tumorType <- tumorType2
      #stop()
    }
    if(sum(refData$metaData$`Tumor type simple` == tumorType) < 5){
      message(paste("Less then 5 reference samples for",tumorType))
    }
  }
  refDataGene <- refData$counts[gene,]
  sampleDataGene <- runData[gene,sample]
  if (!(sample %in% names(refDataGene))){
    dataGene <- c(refDataGene,sampleDataGene)
    names(dataGene)[length(dataGene)] <- sample
    dataGeneOrdered <- dataGene[order(dataGene)]
    metaData <- refData$metaData[names(dataGeneOrdered)[names(dataGeneOrdered) != sample],]
  }else{
    dataGene <- refDataGene
    dataGeneOrdered <- dataGene[order(dataGene)]
    metaData <- refData$metaData[names(dataGeneOrdered),]
  }
  if (showAll){
    #    layout(mat=matrix(ncol=2,nrow=2,data=c(1,2,3,3),byrow = T),heights = c(7,1))
    par(mar=c(5,10,2,10))
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='grey',xlab=paste(gene,"- Counts per million"),ylab="Samples - all tumor types",main=paste(gene,"- All tumor types"),cex=2)
    tumorTypeSamples <- rownames(metaData)[metaData$`Tumor type simple` == tumorType]
    points(dataGeneOrdered[tumorTypeSamples],which(names(dataGeneOrdered) %in% tumorTypeSamples),pch=20,cex=2)
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=3)
    legend("bottomright",legend = c(tumorType,sample),pch=20,col=c('black','red'),bty='n')
  }else{
    par(mar=c(5,10,2,10))
    tumorTypeSamples <- rownames(metaData)[metaData$`Tumor type simple` == tumorType]
    dataGeneOrdered <- dataGeneOrdered[names(dataGeneOrdered) %in% c(tumorTypeSamples,sample)]
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='black',xlab=paste(gene,"- Counts per million"),ylab=paste("Samples -",tumorType),main=paste(gene,"-",tumorType),cex=2)
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=3)
    legend("bottomright",legend = sample,pch=20,col='red',bty='n')
  }
  if (pdf){
    pdf(paste0(baseDirWTS,folder,"/",sample,"_",gene,"_",tumorType,".pdf"), width = 12, height = 7, bg="white")
    layout(mat=matrix(ncol=2,nrow=2,data=c(1,2,3,3),byrow = T),heights = c(7,1))
    par(mar=c(5,4,2,0.5))
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='grey',xlab=paste(gene,"- Counts per million"),ylab="Samples - all tumor types",main=paste(gene,"- All tumor types"),cex=2)
    tumorTypeSamples <- rownames(metaData)[metaData$`Tumor type simple` == tumorType]
    points(dataGeneOrdered[tumorTypeSamples],which(names(dataGeneOrdered) %in% tumorTypeSamples),pch=20,cex=2)
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=3)
    legend("bottomright",legend = c(tumorType,sample),pch=20,col=c('black','red'),bty='n')
    
    dataGeneOrdered <- dataGeneOrdered[names(dataGeneOrdered) %in% c(tumorTypeSamples,sample)]
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='black',xlab=paste(gene,"- Counts per million"),ylab=paste("Samples -",tumorType),main=paste(gene,"-",tumorType),cex=2)
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=3)
    legend("bottomright",legend = sample,pch=20,col='red',bty='n')
    par(mar=c(0,4,0,0))
    plot(0,0,cex=0,xlim=c(0,1),ylim=c(0,1),axes=F,ylab="",xlab="")
    text(x = 0.5,y = 0.5,labels = paste("Reference cohort version:",refData$version))
    dev.off()
  }
  
}



loadRefData <- function(countSet = "20200623_PMCdiag_RNAseq_counts_60357.csv"){
  baseDir <- paste0(baseDirWTS,"QualityControl/expressionData/")
  refFiles <- list.files(baseDir,pattern = "refData.rds")
  
  #  refFileDate <- max(sapply(refFiles, function(x) strsplit(x,"_")[[1]][1]))
  #  refDataDate <- strsplit(countSet,"_")[[1]][1]
  #  if (as.numeric(refFileDate) > as.numeric(refDataDate)){
  if (file.exists(paste0(baseDir,sub(".csv","_refData.rds",countSet)))){
    refData <- readRDS(paste0(baseDir,sub(".csv","_refData.rds",countSet)))
    refData$counts <- apply(refData$counts,2,function(x) (x/sum(x))*1000000)
    return(refData)
  }else{
    countData <- read.csv(paste0(baseDir,countSet),sep="\t",stringsAsFactors = F)
    #countData <- countData[countData$GeneName != "MIR6867",]
    samples <- sapply(colnames(countData)[c(3:ncol(countData))],function(x) strsplit(x,"_")[[1]][1])
    #dups <- samples[duplicated(samples)]
    #samplesDedup <- samples[!(samples %in% dups)]
    samplesDedup <- samples
    
    metaData <- loadRNAseqOverview(samples=samplesDedup,type="biomaterial")
    dupSamples <- rownames(metaData)[grep("\\.1",rownames(metaData))]
    dupSamplesFirst <- sub("\\.1","",dupSamples)
    metaData <- metaData[!(rownames(metaData) %in% dupSamplesFirst),]
    rownames(metaData) <- sub("\\.1","",rownames(metaData))
    
    countDataDedup <- cbind(countData[,c(1,2)],countData[,names(samplesDedup)[samplesDedup %in% rownames(metaData)]])
    colnames(countDataDedup)[c(3:ncol(countDataDedup))] <- sapply(colnames(countDataDedup)[c(3:ncol(countDataDedup))],function(x) strsplit(x,"_")[[1]][1])
    tumorFusion <- metaData[,c("Tumor type simple","Resultaat RNA seq (relevante)","HIX Nr")]
    tumorFusion <- tumorFusion[colnames(countDataDedup)[c(3:ncol(countDataDedup))],]
    
    dupGenes <- countDataDedup$GeneName[duplicated(countDataDedup$GeneName)]
    countDataDedup2 <- countDataDedup[!(countDataDedup$GeneName %in% dupGenes),]
    tempData <- as.data.frame(matrix(ncol = ncol(countDataDedup),nrow = length(unique(dupGenes))))
    colnames(tempData) <- colnames(countDataDedup2)
    tempData[,2] <- unique(dupGenes)
    for ( i in 1:length(unique(dupGenes))){
      tempData[i,c(3:ncol(tempData))] <- apply(countDataDedup[countDataDedup$GeneName == unique(dupGenes)[i],c(3:ncol(countDataDedup))],2,sum)
    }
    countDataDedup <- rbind(countDataDedup2[,-1],tempData[,-1])
    
    rownames(countDataDedup) <- countDataDedup$GeneName
    countDataDedup <- countDataDedup[,-1]
    date <- gsub("-","",Sys.Date())
    tumorFusion <- tumorFusion[!is.na(tumorFusion$`Tumor type simple`),]
    countDataDedup <- countDataDedup[,rownames(tumorFusion)]
    
    mostRecentQC <- list.files(paste(baseDirWTS,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcData <- read.csv(mostRecentQC,sep="\t",stringsAsFactors = F)
    enoughReads <- qcData$Biomaterial.ID[qcData$uniqueReads.10.6. > 30]
    countDataDedup <- countDataDedup[,colnames(countDataDedup) %in% enoughReads]
    
    ## remove wrong samples ##
    countDataDedup[,c("PMABM000AGL","PMABM000AGN")] <- countDataDedup[,c("PMABM000AGN","PMABM000AGL")]
    countDataDedup <- countDataDedup[,!(colnames(countDataDedup) %in% c("PMABM000BJE","PMABM000BJG","PMABM000BJI","PMABM000BJK","PMABM000BJO","PMABM000BJU","PMABM000BKF"))]
    countDataDedup <- countDataDedup[,rownames(tumorFusion)[tumorFusion$`Tumor type simple` !="No HiX diagnosis"]]
    
    tumorFusion <- tumorFusion[colnames(countDataDedup),]
    saveRDS(list("counts" = countDataDedup,"metaData"=tumorFusion,"version"=countSet),paste0(baseDir,sub(".csv","_refData.rds",countSet)))
    
    countDataDedupNorm <- apply(countDataDedup,2,function(x) (x/mean(x))*1000000)
    return(list("counts" = countDataDedupNorm,"metaData"=tumorFusion,"version"=countSet))
  }
}

checkExpressionData <- function(folder){
  expressionFile <- paste0(baseDirWTS,folder,"/expressionData/",folder,"_counts.csv")
  expressionFiles <- getFileList(folder,rootDir = baseDirWTS,pattern = "RNA-Seq.gene_id.exon.counts.txt")$targetFiles
  expressionFiles2 <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",expressionFiles)
  if(!file.exists(expressionFile)){
    message("Expression data not available yet, start downloading...")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }
  if (sum(expressionFiles2 %in% list.files(paste0(baseDirWTS,folder,"/expressionData/"),pattern = "_RNA-Seq.gene_id.exon.counts",full.names = F)) < length(expressionFiles2)){
    message("Expression data not available yet, start downloading...")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }
}

loadRunExpressionData <- function(folder){
  expressionFile <- paste0(baseDirWTS,folder,"/expressionData/",folder,"_counts.csv")
  expressionFiles <- getFileList(folder,rootDir = baseDirWTS,pattern = "RNA-Seq.gene_id.exon.counts.txt")$targetFiles
  expressionFiles2 <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",expressionFiles)
  if(!file.exists(expressionFile)){
    message("Expression data not available yet, start downloading...")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }
  if (sum(expressionFiles2 %in% list.files(paste0(baseDirWTS,folder,"/expressionData/"),pattern = "_RNA-Seq.gene_id.exon.counts",full.names = F)) < length(expressionFiles2)){
    message("Expression data not available yet, start downloading...")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }else{
    counts <- read.csv(expressionFile,sep="\t",stringsAsFactors = F)
  }
  return(apply(counts,2,function(x) (x/sum(x))*1000000))
}

downloadExpressionData <- function(folder){
  expressionFiles <- getFileList(folder,rootDir = baseDirWTS,pattern = "RNA-Seq.gene_id.exon.counts.txt")$targetFiles
  dir.create(paste0(baseDirWTS,folder,"/expressionData"),showWarnings = F)
  
  dataList <- list()
  for ( i in 1:length(expressionFiles)){
    expFileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",expressionFiles[i])
    destFile <- paste(baseDirWTS,folder,"/expressionData/",expFileName,sep="")
    if ( file.exists(destFile)){
      message(paste(strsplit(expFileName,"_")[[1]][1] ,"already downloaded"))
    }else{
      message(paste("downloading",expFileName))
      GET(expressionFiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
    }
    dataList[[i]] <- read.csv(destFile,sep="\t",stringsAsFactors = F,skip=3)
  }
  countData <- as.data.frame(do.call(cbind,lapply(dataList,function(x) as.numeric(x$Counts))))
  countData$GeneName <- dataList[[1]]$GeneName
  colnames(countData) <- c(sapply(expressionFiles,function(x) strsplit(gsub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",x),"_")[[1]][1] ),"GeneName")
  dupGenes <- countData$GeneName[duplicated(countData$GeneName)]
  countData2 <- countData[!(countData$GeneName %in% dupGenes),]
  tempData <- as.data.frame(matrix(ncol = ncol(countData),nrow = length(unique(dupGenes))))
  colnames(tempData) <- colnames(countData2)
  tempData[,ncol(tempData)] <- unique(dupGenes)
  for ( i in 1:length(unique(dupGenes))){
    tempData[i,c(1:(ncol(countData)-1))] <- apply(countData[countData$GeneName == unique(dupGenes)[i],c(1:(ncol(countData)-1))],2,sum)
  }
  countData <- rbind(countData2,tempData)
  
  rownames(countData) <- countData$GeneName
  countData <- countData[,-ncol(countData)]
  return(countData)
}

