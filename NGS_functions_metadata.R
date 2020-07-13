makeMetaDataWES <- function(seqRunDir=NULL,rootDir=baseDirWES){
  date <- Sys.Date()
  date <- gsub("-","_",date)
  
  dataLinks <- getFileList(seqRunDir,rootDir,"WXS.qci.vcf")
  if(identical(dataLinks,1)){
    return(1)
  }else{
    vcfFiles <- dataLinks$targetFiles
    for ( i in 1:length(vcfFiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",vcfFiles[i])
      destFile <- paste(rootDir,seqRunDir,"/",fileName,sep="")
      GET(vcfFiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
    }
    
    multiQCfiles <- getFileList(seqRunDir, rootDir, "multiqc_data.zip")$targetFiles
    
    downloadedSamples <- list.files(paste0(rootDir,"QualityControl/multiQCfiles/"),pattern = "zip",full.names = T)
    
    for ( i in 1:length(multiQCfiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",multiQCfiles[i])
      destFile <- paste(rootDir,"QualityControl/multiQCfiles/",fileName,sep="")
      if ( !(destFile %in% downloadedSamples)){
        GET(multiQCfiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
      }
      destDir <- sub("_WXS.multiqc_data.zip","",destFile)
      if ( !dir.exists(destDir)){
        dir.create(destDir,showWarnings = F)  
        unzip(destFile,exdir = destDir)
      }else if(!dir.exists(paste(destDir,"/",sub(".zip","",fileName),sep=""))){
        unzip(destFile,exdir = destDir)
      }
    }
    samples <- sub("_WXS.multiqc_data.zip","",multiQCfiles)
    samples <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",samples)
    
    samples <- samples[grep("zip",samples,invert = T)]
    
    pairedSamples <- samples[sapply(samples,function(x) length(strsplit(x,"_")[[1]])) == 3]
    singleSamples <- samples[sapply(samples,function(x) length(strsplit(x,"_")[[1]])) == 2]
    
    dataTypesSingle <- c("fastqc","gatk_varianteval","general_stats","picard_AlignmentSummaryMetrics","picard_baseContent","picard_dups","picard_gcbias","picard_insertSize","picard_HsMetrics","picard_varientCalling","verifybamid")
    dataTypesPaired <- c("gatk_varianteval","general_stats","picard_OxoGMetrics","tumor-normal_comparison")
    dataTypesSingle2 <- paste0("multiqc_",dataTypesSingle,".txt")
    dataTypesPaired2 <- paste0("multiqc_",dataTypesPaired,".txt")
    
    ## load data from multiQC files ##
    
    dataList <- list()
    for ( i in 1:length(pairedSamples)){
      curSampleDir <- list.dirs(paste0(rootDir,"QualityControl/multiQCfiles/",pairedSamples[i]),recursive = T)[2]
      if (is.na(curSampleDir)){
        next
      }
      tempList <- list()
      for ( j in 1:length(dataTypesPaired2)){
        if(file.exists(paste0(curSampleDir,"/",dataTypesPaired2[j]))){
          tempList[[j]] <- read.csv(paste0(curSampleDir,"/",dataTypesPaired2[j]),sep="\t",stringsAsFactors = F)
          colnames(tempList[[j]]) <- paste(dataTypesPaired[j],colnames(tempList[[j]]),sep="_")
        }
      }
      dataList[[i]] <- as.data.frame(t(unlist(tempList)))
    }
    
    colnamesAll <- unique(unlist(lapply(dataList, function(x) colnames(x))))
    
    multiQCdatapaired <- as.data.frame(matrix(nrow=length(dataList),ncol=length(colnamesAll)))
    colnames(multiQCdatapaired) <- colnamesAll
    for ( i in 1:length(pairedSamples)){
      if ( !is.null(dataList[[i]])){
        multiQCdatapaired[i,colnames(dataList[[i]])] <- as.matrix(dataList[[i]])   
      }
    }
    rownames(multiQCdatapaired) <- pairedSamples
    multiQCdatapaired[,"Biomaterial ID pair"] <- sapply(rownames(multiQCdatapaired), function(x) paste(strsplit(x,"_")[[1]][c(1,2)],collapse="_"))
    
    dups <- multiQCdatapaired$`Biomaterial ID`[duplicated(multiQCdatapaired$`Biomaterial ID pair`)]
    multiQCdatapaired <- multiQCdatapaired[!(multiQCdatapaired$`Biomaterial ID pair` %in% dups),]
    multiQCdatapaired$sample1 <- sapply(multiQCdatapaired$`Biomaterial ID pair`,function(x) strsplit(x,"_")[[1]][1])
    multiQCdatapaired$sample2 <- sapply(multiQCdatapaired$`Biomaterial ID pair`,function(x) strsplit(x,"_")[[1]][2])
    multiQCdatapaired <- multiQCdatapaired[,grep("OxoG",colnames(multiQCdatapaired),invert = T)]
    colnames(multiQCdatapaired) <- paste("pair_",colnames(multiQCdatapaired),sep="")
    
    singleSamplesBM <- unique(c(multiQCdatapaired$pair_sample1,multiQCdatapaired$pair_sample2))
    missingSamples <- singleSamplesBM[!(singleSamplesBM %in% sapply(singleSamples,function(x) strsplit(x,"_")[[1]][1]))]
    if(length(missingSamples) > 0){
      sampleDirs <- list.dirs(paste0(rootDir,"QualityControl/multiQCfiles/"),full.names = F)
      singleSamples <- c(singleSamples,sampleDirs[grep(paste0("^",missingSamples[1],"_PMCRZ"),sampleDirs)][1])
    }
    singleSamples <- singleSamples[sapply(singleSamples,function(x) strsplit(x,"_")[[1]][1]) %in% singleSamplesBM]
    dataList <- list()
    for ( i in 1:length(singleSamples)){
      curSampleDir <- list.dirs(paste0(rootDir,"QualityControl/multiQCfiles/",singleSamples[i]),recursive = T)[2]
      if (is.na(curSampleDir)){
        next
      }
      tempList <- list()
      for ( j in 1:length(dataTypesSingle2)){
        tempList[[j]] <- read.csv(paste0(curSampleDir,"/",dataTypesSingle2[j]),sep="\t",stringsAsFactors = F)
        colnames(tempList[[j]]) <- paste(dataTypesSingle[j],colnames(tempList[[j]]),sep="_")
      }
      dataList[[i]] <- as.data.frame(t(unlist(tempList)))
    }
    
    colnamesAll <- unique(unlist(lapply(dataList, function(x) colnames(x))))
    
    multiQCdataSingle <- as.data.frame(matrix(nrow=length(dataList),ncol=length(colnamesAll)))
    colnames(multiQCdataSingle) <- colnamesAll
    for ( i in 1:length(singleSamples)){
      if ( !is.null(dataList[[i]])){
        multiQCdataSingle[i,colnames(dataList[[i]])] <- as.matrix(dataList[[i]])   
      }
    }
    rownames(multiQCdataSingle) <- sapply(singleSamples,function(x) paste(strsplit(x,"_")[[1]][c(1,2)],collapse="_"))
    multiQCdataSingle[,"Biomaterial ID"] <- sapply(multiQCdataSingle$fastqc_Sample, function(x) strsplit(x,"_")[[1]][1])
    
    dups <- multiQCdataSingle$`Biomaterial ID`[duplicated(multiQCdataSingle$`Biomaterial ID`)]
    multiQCdataSingle <- multiQCdataSingle[!(multiQCdataSingle$`Biomaterial ID` %in% dups),]
    
    tempList <- list()
    for ( i in 1:nrow(multiQCdataSingle)){
      tempList[[i]] <- as.vector(c(multiQCdataSingle[i,],multiQCdatapaired[grep(multiQCdataSingle$`Biomaterial ID`[i],multiQCdatapaired$`pair_Biomaterial ID pair`)[1],]))
    }
    cramFiles <- getFileList(seqRunDir,rootDir,pattern = ".cram")$targetFiles
    names(cramFiles) <- sapply(sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",cramFiles), function(x) strsplit(x,"_")[[1]][1])
    multiQCdata <- do.call(rbind,tempList)
    multiQCdata <- as.data.frame(apply(multiQCdata,2,unlist),stringsAsFactors = F)
    rownames(multiQCdata) <- multiQCdata$`Biomaterial ID`
    samples <- c(multiQCdatapaired$pair_sample1,multiQCdatapaired$pair_sample2)
    if (!sum(samples %in% multiQCdata$`Biomaterial ID`) == length(samples)){
      missingSamples <- samples[!samples %in% multiQCdata$`Biomaterial ID`]
      missingSampleWESoverview <- loadWESOverview(samples = missingSamples)
      mostRecentQC <- list.files(paste(rootDir,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
      mostRecentQC <- mostRecentQC[length(mostRecentQC)]
      qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
      for( i in 1:length(missingSamples)){
        multiQCdata <- rbind(multiQCdata,qcdataAll[qcdataAll$`Biomaterial ID` == missingSamples[i],colnames(multiQCdata)])
        missingCramFiles <- getFileList(seqRunDir = missingSampleWESoverview$seqRunID[missingSampleWESoverview$`Biomaterial ID` == missingSamples[i]][1],rootDir = rootDir,pattern = ".cram")$targetFiles
        names(missingCramFiles) <- sapply(sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",missingCramFiles), function(x) strsplit(x,"_")[[1]][1])
        cramFiles <- c(cramFiles,missingCramFiles[missingSamples[i]])
      }
      
      
    }
    
    #  BMlist <- as.data.frame(as.matrix(read.xlsx(fileList$bmlijst)),stringsAsFactors = F)
    #  BMinfo <- BMlist[BMlist$PMABM %in% samples,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer","DIN")]
    NGSinfo <- as.data.frame(as.matrix(read.xlsx(fileList$NGSdiagnostiek,sheet = "DNA library prep")),stringsAsFactors = F)
    NGSinfo <- NGSinfo[NGSinfo$PMABM %in% samples,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer","Gem.Fragment.lengte.BA","nM.DNA","totaal.in.oprep","DIN","totaal.(ng)")]
    dups <- rev(duplicated(rev(NGSinfo$PMABM)))
    NGSinfo <- NGSinfo[!dups,]
    otherSamples <- samples[!(samples %in% NGSinfo$PMABM)]
    otherMatrix <- as.data.frame(matrix(ncol=ncol(NGSinfo),nrow=length(otherSamples)))
    colnames(otherMatrix) <- colnames(NGSinfo)
    otherMatrix$PMABM <- otherSamples
    otherMatrix$PMABS <- otherSamples
    NGSinfo <- rbind(NGSinfo,otherMatrix)
    rownames(NGSinfo) <- NGSinfo$PMABM
    multiQCdata <- cbind(multiQCdata,NGSinfo[multiQCdata$`Biomaterial ID`,])
    multiQCdata$insertSize_BA <- as.numeric(multiQCdata$Gem.Fragment.lengte.BA)
    multiQCdata$concentration_nM <- round(as.numeric(multiQCdata$nM.DNA),2)
    multiQCdata$yield <- round(as.numeric(multiQCdata$`totaal.(ng)`,0))
    
    
    
    # NGSlist <- read.xlsx(fileList$NGSdiagnostiek,sheet = "DNA library prep")
    # libData <- NGSlist[,c("PMABM","Gem.Fragment.lengte.BA","nM.DNA","totaal.in.oprep","DIN","totaal.(ng)")]
    # colnames(libData) <- c("Biomaterial ID","insertSize_BA","concentration_nM","input_ng","DIN","yield(ng)")
    # libData$DIN <- round(as.numeric(libData$DIN),2)
    # libData <- libData[!is.na(libData$`Biomaterial ID`),]
    # libData <- libData[libData$`Biomaterial ID` %in% samples,]
    # dups <- rev(duplicated(rev(libData$`Biomaterial ID`)))
    # libData <- libData[!dups,]
    # rownames(libData) <- libData$`Biomaterial ID`
    # multiQCdata <- cbind(multiQCdata,libData[multiQCdata$`Biomaterial ID`,])
    # 
    # 
    # BSlijst <- read.xlsx(fileList$bslijst)
    # BSlijst <- BSlijst[BSlijst$PMABS %in% multiQCdata$PMABS,c("Tumor.%","PMABS"),]
    # BSlijst <- BSlijst[!is.na(BSlijst$PMABS),]
    # otherBS <- multiQCdata$PMABS[!(multiQCdata$PMABS %in% BSlijst$PMABS)]
    # if (length(otherBS) > 0){
    #   otherBSMat <- matrix(ncol=2,nrow=length(otherBS))
    #   colnames(otherBSMat) <- colnames(BSlijst)
    #   otherBSMat[,"PMABS"] <- otherBS
    #   BSlijst <- rbind(BSlijst,otherBSMat)
    # }
    # rownames(BSlijst) <- BSlijst$PMABS
    
    # multiQCdata <- cbind(multiQCdata,BSlijst[multiQCdata$PMABS,"Tumor.%"])
    # colnames(multiQCdata)[ncol(multiQCdata)] <- "tumorPerc"
    # multiQCdata$tumorPerc <- as.character(multiQCdata$tumorPerc)
    multiQCdata$tumorPerc <- rep(NA,nrow(multiQCdata))
    multiQCdata$uniqueReads <- round(as.numeric(multiQCdata$fastqc_Total.Sequences) * as.numeric(multiQCdata$fastqc_total_deduplicated_percentage) / 100000000,0)
    itherList <- as.data.frame(read.xlsx(fileList$ither),stringsAsFactors=F)
    itherList <- itherList[!is.na(itherList$HIX),]
    itherList$PMABS.tumor[is.na(itherList$PMABS.tumor)] <- "unknown"
    metaData <- as.data.frame(matrix(nrow=nrow(multiQCdatapaired),ncol=19))
    colnames(metaData) <- c("Skion ID","HiX Nr","PMABS tumor","PMABM tumor","PMABS normal","PMABM normal","T-nummer","Vraagstelling","Tumor perc","UniqueReadsTumor(10^6)","UniqueReadsNormal(10^6)","MeanCoverageTumor","MeanCoverageNormal","novelVariants","ContaminationTumor","ContaminationNormal","TumorNormalCheck","cramTumor","cramNormal")
    for ( i in 1:nrow(metaData)){
      metaData[i,c(1,2,3,4,7,9,10,12,15,17)] <- multiQCdata[multiQCdata$`Biomaterial ID` == multiQCdatapaired$pair_sample1[i],c("PMCBS","HIX","PMABS","PMABM","PA-nummer","tumorPerc","uniqueReads","picard_HsMetrics_MEAN_TARGET_COVERAGE","verifybamid_FREEMIX","pair_tumor-normal_comparison_Ratio.as.expected.")]
      metaData[i,c(5,6,11,13,16)] <- multiQCdata[multiQCdata$`Biomaterial ID` == multiQCdatapaired$pair_sample2[i],c("PMABS","PMABM","uniqueReads","picard_HsMetrics_MEAN_TARGET_COVERAGE","verifybamid_FREEMIX")]
      metaData[i,14] <- multiQCdatapaired[i,"pair_gatk_varianteval_novel_sites"]
      metaData[i,c(18,19)] <- cramFiles[c(multiQCdatapaired$pair_sample1[i],multiQCdatapaired$pair_sample2[i])]
      if(metaData[i,"PMABS tumor"] %in% itherList$PMABS.tumor){
        metaData[i,"Vraagstelling"] <- paste("iTHER",itherList$Ither.nummer[itherList$PMABS.tumor == metaData[i,"PMABS tumor"]])
      }else{
        metaData[i,"Vraagstelling"] <- "Diagnostiek"
      }
    }
    metaData$Vraagstelling[grep("PMRBM",metaData[,"PMABM tumor"])] <- "Research"
    metaData$Vraagstelling[grep("PMOBM",metaData[,"PMABM tumor"])] <- "Organoid"
    metaData$MeanCoverageTumor <- round(as.numeric(metaData$MeanCoverageTumor),0)
    metaData$MeanCoverageNormal <- round(as.numeric(metaData$MeanCoverageNormal),0)
    write.xlsx(metaData,paste0(rootDir,seqRunDir,"/metaData_LKR.xlsx"))
    
    vcfFiles2 <- getFileList(seqRunDir,rootDir,pattern = "WXS.vcf.gz")$targetFiles
    if(!dir.exists(paste0(rootDir,seqRunDir,"/rawData"))){
      dir.create(paste0(rootDir,seqRunDir,"/rawData"))
    }
    for ( i in 1:length(vcfFiles2)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",vcfFiles2[i])
      destFile <- paste(rootDir,seqRunDir,"/rawData/",fileName,sep="")
      if (!file.exists(destFile)){
        GET(vcfFiles2[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = F))
        if(length(strsplit(fileName,"_")[[1]]) == 3){
          alissaSingleSampleVCF(folder = seqRunDir, vcfFile = fileName)
        }else if(length(strsplit(fileName,"_")[[1]]) == 4){
          alissaSomaticVCF(folder = seqRunDir, vcfFile = fileName)
        }
      }
    }
    
    
    return(metaData)
  }
}






makeFusionExcelFiles <- function(seqRunDir,rootDir = baseDirWTS){
  date <- Sys.Date()
  date <- gsub("-","_",date)
  
  dataLinks <- getFileList(seqRunDir,rootDir,"star-fusion_predicted.annotated.filtered.tsv")
  if(identical(dataLinks,1)){
    return(1)
    
  }else{
    fusionFiles <- dataLinks$targetFiles
    samples <- sapply(sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",fusionFiles),function(x) strsplit(x,"_")[[1]][1])
    
    dataList <- list()
    multiQCfiles <- getFileList(seqRunDir,rootDir,"multiqc_data.zip")$targetFiles
    for ( i in 1:length(multiQCfiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",multiQCfiles[i])
      destFile <- paste(baseDirWTS,"QualityControl/multiQCfiles/",fileName,sep="")
      GET(multiQCfiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
      destDir <- sub("_RNA-Seq.multiqc_data.zip","",destFile)
      if ( !dir.exists(destDir)){
        dir.create(destDir,showWarnings = F)  
        unzip(destFile,exdir = destDir)
      }else if(!dir.exists(paste(destDir,"/",fileName,sep=""))){
        unzip(destFile,exdir = destDir)
      }  
      
      dataList[[i]] <- as.data.frame((read.csv(paste0(list.dirs(destDir)[2],"/","multiqc_fastqc.txt"),sep="\t",stringsAsFactors = F)))
    }
    
    maxColSample <- which.max(unlist(lapply(dataList,function(x) dim(x)[2])))
    colnamesAll <- unique(unlist(lapply(dataList, function(x) colnames(x))))
    
    multiQCdata <- as.data.frame(matrix(nrow=length(dataList),ncol=length(colnamesAll)))
    colnames(multiQCdata) <- colnamesAll
    for ( i in 1:length(samples)){
      if ( !is.null(dataList[[i]])){
        multiQCdata[i,colnames(dataList[[i]])] <- as.matrix(dataList[[i]])   
      }
    }
    multiQCdata$PMABM <- sapply(multiQCdata$Sample,function(x) strsplit(x,"_")[[1]][1])
    rownames(multiQCdata) <- multiQCdata$PMABM
    
    
    NGSlist <- suppressWarnings(read.xlsx(fileList$NGSdiagnostiek,sheet = "RNA Library prep"))
    NGSlist$Opmerking[is.na(NGSlist$Opmerking)] <- "geenopmering"
    NGSlist <- NGSlist[NGSlist$Opmerking != "UDI NextFlex",]
    libData <- NGSlist[,c("PMCBS","PMABS","PMABM","HIX","Type.BS","PA-nummer","Gem.Fragment.lengte.BA","nM.DNA","Volume.(ul)","totaal.in.oprep","RIN")]
    colnames(libData) <- c("PMCBS","PMABS","Biomaterial ID","HIX","Type.BS","PA-nummer","insertSize_BA","concentration_nM","volume_ul","input_ng","RIN")
    libData$yield <- as.numeric(libData$concentration_nM) * as.numeric(libData$volume_ul)
    libData$RIN <- round(as.numeric(libData$RIN),2)
    libData <- libData[!is.na(libData$`Biomaterial ID`),]
    libData <- libData[libData$`Biomaterial ID` %in% samples,]
    dups <- rev(duplicated(rev(libData$`Biomaterial ID`)))
    libData <- libData[!dups,]
    
    rownames(libData) <- libData$`Biomaterial ID`
    
    # BMlist <- as.data.frame(as.matrix(read.xlsx(fileList$bmlijst)),stringsAsFactors = F)
    # datum.isolatie <- convertToDate(BMlist[,"Datum.isolatie"])
    # dd <- cbind(as.data.frame(BMlist[,c("PMABS","PMABM","RIN")]),datum.isolatie)
    # colnames(dd)[2] <- "Biomaterial ID"
    # dd <- dd[!is.na(dd$`Biomaterial ID`),]
    # dd <- dd[dd$`Biomaterial ID` %in% samples,]
    # rownames(dd) <- dd$`Biomaterial ID`
    cramFiles <- getFileList(seqRunDir,baseDirWTS,pattern = "RNA-Seq.cram")$targetFiles
    names(cramFiles) <- sapply(sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",cramFiles),function(x) strsplit(x,"_")[[1]][1])
    libData <- cbind(libData,as.numeric(multiQCdata[rownames(libData),"Total.Sequences"]) * as.numeric(multiQCdata[rownames(libData),"total_deduplicated_percentage"])/100,as.numeric(multiQCdata[rownames(libData),"Total.Sequences"]),cramFiles[rownames(libData)])
    colnames(libData)[c((ncol(libData)-2):ncol(libData))] <- c("unique_reads","total_reads","cramFile")
    
    blackList <- read.csv(fileList$fusionBlacklist,sep="\t",stringsAsFactors = F)
    
    samplesDone <- c()
    metaList <- list()
    
    perLine2 <- dataLinks$perLine
    dir.create(paste0(baseDirWTS,seqRunDir,"/tsvFiles/"),showWarnings = F)
    for ( i in 1:length(fusionFiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",fusionFiles[i])
      destFile <- paste0(baseDirWTS,seqRunDir,"/tsvFiles/",fileName)
      GET(fusionFiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
      
      
      sampleName <- strsplit(fileName,"_")[[1]][1]
      topHeader <- perLine2[grep(paste0(sampleName," tumor"),perLine2):(grep(paste0(sampleName," tumor"),perLine2)+3)]
      curFileHeader <- t(t(as.matrix(read.csv(destFile,sep="\t",stringsAsFactors = F,header=F))))
      rowsToSkip <- grep("FusionName",(curFileHeader)[,1]) - 1
      curFileHeader <- as.vector(curFileHeader[c(1:rowsToSkip),1])
      curMetaData <- libData[sampleName,]
      metaList[[i]] <- curMetaData
      libraryMetaData <- matrix(ncol=1,nrow=6)
      libraryMetaData[1,1] <- paste("## PMABS: ",curMetaData$PMABS,sep="")
      libraryMetaData[2,1] <- paste("## RIN: ",curMetaData$RIN,sep="")
      libraryMetaData[3,1] <- paste("## Input ng: ",curMetaData$input_ng,sep="")
      libraryMetaData[4,1] <- paste("## Insert size BA: ",curMetaData$insertSize_BA,sep="")
      libraryMetaData[5,1] <- paste("## Yield (fmol): ",round(curMetaData$yield,0),sep="")
      libraryMetaData[6,1] <- paste("## Unique reads: ",round(curMetaData$unique_reads,0),sep="")
      curFileHeader <- c(curFileHeader,libraryMetaData)
      curFileData <- read.csv(destFile,sep="\t",stringsAsFactors = F,header=F,skip=rowsToSkip)
      blackListTable <- as.data.frame(matrix(ncol=3,nrow=nrow(curFileData)),stringsAsFactors = F)
      blackListTable[1,] <- c("timesObserved","relevant","note")
      if (nrow(curFileData) > 1){
        for ( j in 2:nrow(curFileData)){
          for ( k in 1:ncol(curFileData)){
            if (!is.na(curFileData[j,k]) & nchar(curFileData[j,k]) > 32766){
              curFileData[j,k] <- "Too large for excel to display, check tsv file"
            }
          }
        }
        for ( j in 2:nrow(curFileData)){
          genes <- strsplit(curFileData[j,1],"--")[[1]]
          if (sum(genes[1] == blackList$gene1 & genes[2] == blackList$gene2 ) > 0){
            blackListTable[j,] <- blackList[which(genes[1] == blackList$gene1 & genes[2] == blackList$gene2),c("observed","Relevant","note")]
          }
          else if (sum(genes[2] == blackList$gene1 & genes[1] == blackList$gene2 ) > 0){
            blackListTable[j,] <- blackList[which(genes[2] == blackList$gene1 & genes[1] == blackList$gene2),c("observed","Relevant","note")]
          }
          else{
            blackListTable[j,] < c(0,NA,NA)
          }
        }
      }
      curFileData <- cbind(curFileData[,c(1:10)],blackListTable,curFileData[,c(11:ncol(curFileData))])
      topHeader <- c("Resultaten RNAseq voor Fusiegendetectie"," ",topHeader)
      totalRows <- length(topHeader) + length(curFileHeader) + nrow(curFileData) + 1
      totalCols <- ncol(curFileData)
      totalMatrix <- matrix(nrow=totalRows,ncol = totalCols)
      totalMatrix[c(1:length(topHeader)),1] <- topHeader
      totalMatrix[c((length(topHeader)+2):(length(topHeader)+1+length(curFileHeader))),1] <- curFileHeader
      totalMatrix[c((length(topHeader)+2+length(curFileHeader)):nrow(totalMatrix)),] <- as.matrix(curFileData)
      if(identical(grep(sampleName,list.files(paste0(baseDirWTS,seqRunDir),pattern = ".xlsx")),integer(0))){
        tryCatch(write.xlsx(totalMatrix,paste0(baseDirWTS,seqRunDir,"/",sampleName,"_LKR.xlsx"),colNames = F,overwrite=F),error = function(e) NULL ) 
      }
      samplesDone <- c(samplesDone,sampleName)
      #print(paste("Downloaded sample:",sampleName,"and generated excel file"))
    }
    
    #print("Done!")
    
    samples <- libData[libData$`Biomaterial ID` %in% samplesDone,c("PMCBS","HIX","PMABS","Biomaterial ID","Type.BS","PA-nummer","RIN")]
    
    #    BSlijst <- read.xlsx(fileList$bslijst)
    #    BSlijst <- BSlijst[BSlijst$PMABS %in% samples$PMABS,c("Tumor.%","PMABS")]
    #    rownames(BSlijst) <- make.unique(BSlijst$PMABS)
    itherList <- as.data.frame(read.xlsx(fileList$ither),stringsAsFactors=F)
    itherList <- itherList[!is.na(itherList$HIX),]
    itherList$PMABS.tumor[is.na(itherList$PMABS.tumor)] <- "unknown"
    vraagStelling <- rep("Diagnostiek",nrow(samples))
    for( i in 1:nrow(samples)){
      if(samples$PMABS[i] %in% itherList$PMABS.tumor){
        vraagStelling[i] <- paste("iTHER",itherList$Ither.nummer[itherList$PMABS.tumor == samples$PMABS[i]])
      }
    }
    #    samples <- as.data.frame(cbind(samples[,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer")],vraagStelling,BSlijst[samples$PMABS,"Tumor.%"],samples[,c("RIN")]),stringsAsFactors=F)
    samples <- as.data.frame(cbind(samples[,c("PMCBS","HIX","PMABS","Biomaterial ID","Type.BS","PA-nummer")],vraagStelling,rep(NA,nrow(samples)),samples[,c("RIN")]),stringsAsFactors=F)
    
    metaData <- do.call(rbind,metaList)
    samples$yield <- round(metaData[samples$`Biomaterial ID`,"yield"],0)
    samples$uniqueReads <- round(metaData[samples$`Biomaterial ID`,"unique_reads"]/1000000,0)
    samples$input <- metaData[samples$`Biomaterial ID`,"input_ng"]
    samples$cramFiles <- as.character(metaData[samples$`Biomaterial ID`,"cramFile"])
    colnames(samples) <- c("SKION ID","HIX Nr","BioSource ID","Biomaterial ID","Materiaal","Weefsel#","Vraagstelling","TumorPercentage","RIN","yield(fmol)","uniqueReads(10^6)","input(ng)","cramFile")
    #samples$RIN <- round(as.numeric(as.character(samples$RIN)),1)
    
    
    otherSamples <- samplesDone[!(samplesDone %in% libData$`Biomaterial ID`)]
    if(length(otherSamples) > 0){
      otherMatrix <- as.data.frame(matrix(ncol=ncol(samples),nrow=length(otherSamples)))
      colnames(otherMatrix) <- colnames(samples)
      otherMatrix$`Biomaterial ID` <- otherSamples
      otherMatrix$Vraagstelling[grep("PMRBM",otherMatrix$`Biomaterial ID`)] <- "Research"
      otherMatrix$Vraagstelling[grep("PMOBM",otherMatrix$`Biomaterial ID`)] <- "Organoid"
      NGSlist <- NGSlist[!is.na(NGSlist$PMABM),]
      for ( i in 1:nrow(otherMatrix)){
        otherMatrix[i,c("SKION ID","HIX Nr","BioSource ID","Materiaal","Weefsel#","RIN","input(ng)")] <- NGSlist[NGSlist$PMABM == otherMatrix$`Biomaterial ID`[i],c("PMCBS","HIX","PMABS","Type.BS","PA-nummer","RIN","totaal.in.oprep")]
        otherMatrix[i,c("yield(fmol)","uniqueReads(10^6)","cramFile")] <- c(round(metaData[otherMatrix$`Biomaterial ID`[i],"yield"],0),round(metaData[otherMatrix$`Biomaterial ID`[i],"unique_reads"]/1000000,0),cramFiles[otherMatrix$`Biomaterial ID`[i]])
      }
      samples <- rbind(samples,otherMatrix)
    }
    samples$RIN <- as.numeric(gsub("[^0-9\\.]", "", as.character(samples$RIN)))
    #samples$TumorPercentage <- as.numeric(gsub("[^0-9\\.]", "", as.character(samples$TumorPercentage)))
    
    tryCatch(write.xlsx(samples,paste0(baseDirWTS,seqRunDir,"/metaData_LKR.xlsx"),colNames = T,overwrite=T),error = function(e) print("metadata file already exists") ) 
    #return(paste("Created excelsheets and metadata for",nrow(samples),"samples"))
    return(samples)
  }
}

alissaSingleSampleVCF <- function(folder,vcfFile,addAF=F){
  fullVcfFile <- paste0(baseDirWES,folder,"/rawData/",vcfFile)
  vcfAlissa <- sub(".vcf.gz","_Alissa.vcf",fullVcfFile)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  gunzip(fullVcfFile)
  fullVcfFile <- sub(".gz","",fullVcfFile)
  vcf <- vcfR::read.vcfR(fullVcfFile)
  
  vcf@meta <- c("##reference=GRCh38",vcf@meta)
  
  altAlleles <- sapply(vcf@fix[,"ALT"],function(x) strsplit(x,",")[[1]])
  moreThan4 <- which(lapply(altAlleles,length) > 4)
  if(length(moreThan4)>0){
    for ( i in 1:length(moreThan4)){
      fixed <- vcf@fix[moreThan4[i],]
      gt <- strsplit(vcf@gt[moreThan4[i],3],":")[[1]]
      gtAD <- unlist(strsplit(gt[2],",")[[1]])
      gtGT <- unlist(strsplit(gt[1],"/")[[1]])
      fixedAlt <- unlist(strsplit(fixed["ALT"],","))
      allelesToKeep <- (order(as.numeric(gtAD)[-1],decreasing = T)[c(1:4)])
      allelesToKeep <- allelesToKeep[order(allelesToKeep)]
      fixed["ALT"] <- fixedAlt[allelesToKeep]
      gtAD <- gtAD[c(1,allelesToKeep+1)]
      paste(fixedAlt,collapse=",")
      fixed["FILTER"] <- paste0(fixed["FILTER"],";more_than_4_alt_alleles")
      gtGT[2] <- "4"
      gt[1] <- paste(gtGT,collapse = "/")
      gt[2] <- paste(gtAD,collapse = ",")
      vcf@gt[moreThan4[i],3] <- paste(gt,collapse = ":")
      vcf@fix[moreThan4[i,]] <- fixed
    }
  }
  if(addAF){
    gtOptions <- gtOptionsAF <- unique(vcf@gt[,1])
    for (i in 1:length(gtOptions)){
      FORMAT <- strsplit(gtOptions[i],":AD:")[[1]]
      gtOptionsAF[i] <- paste0(FORMAT[1],":AD:AF:",FORMAT[2])
    }
    names(gtOptionsAF) <- gtOptions
    vcf@gt[,1] <- gtOptionsAF[vcf@gt[,1]]
    gt <- t(apply(vcf@gt,1,function(x) {
      #x[1] <- gtOptionsAF[x[1]]
      if(!is.na(x[2])){
        gtSNV <- strsplit(x[2],":")[[1]]
        gtSNVaf <- as.numeric(strsplit(gtSNV[2],",")[[1]])
        gtSNVaf <- paste(round(gtSNVaf[-1]/sum(gtSNVaf),2),collapse=",")
        gtSNV <- c(gtSNV[c(1,2)],gtSNVaf,gtSNV[c(3:length(gtSNV))])
        gtSNV <- paste(gtSNV,collapse = ":")
        x[2] <- gtSNV
      }else{
        gtSNV <- strsplit(x[3],":")[[1]]
        gtSNVaf <- as.numeric(strsplit(gtSNV[2],",")[[1]])
        gtSNVaf <- paste(round(gtSNVaf[-1]/sum(gtSNVaf),2),collapse=",")
        gtSNV <- c(gtSNV[c(1,2)],gtSNVaf,gtSNV[c(3:length(gtSNV))])
        gtSNV <- paste(gtSNV,collapse = ":")
        x[3] <- gtSNV
      }
      return(x)
    }))
    vcf@gt <- gt
  }
  vcfR::write.vcf(vcf,vcfAlissa)
  gzip(fullVcfFile)
  
}

alissaSomaticVCF <- function(folder,vcfFile){
  fullVcfFile <- paste0(baseDirWES,folder,"/rawData/",vcfFile)
  vcfAlissa <- sub(".vcf.gz","_Alissa.vcf",fullVcfFile)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  gunzip(fullVcfFile)
  fullVcfFile <- sub(".gz","",fullVcfFile)
  vcf <- vcfR::read.vcfR(fullVcfFile)
  tumor <- strsplit(vcfFile,"_")[[1]][1]
  normal <- strsplit(vcfFile,"_")[[1]][2]
  
  vcf@meta <- c(vcf@meta[1],"##reference=GRCh38",vcf@meta[2:length(vcf@meta)])
  
  altAlleles <- sapply(vcf@fix[,"ALT"],function(x) strsplit(x,",")[[1]])
  moreThan4 <- which(lapply(altAlleles,length) > 4)
  for ( i in 1:length(moreThan4)){
    fixed <- vcf@fix[moreThan4[i],]
    gtTumor <- strsplit(vcf@gt[moreThan4[i],tumor],":")[[1]]
    gtNormal <- strsplit(vcf@gt[moreThan4[i],normal],":")[[1]]
    gtTumorAF <- unlist(strsplit(gtTumor[3],",")[[1]])
    gtNormalAF <- unlist(strsplit(gtNormal[3],",")[[1]])
    gtTumorAD <- unlist(strsplit(gtTumor[2],",")[[1]])
    gtNormalAD <- unlist(strsplit(gtNormal[2],",")[[1]])
    gtTumorGT <- unlist(strsplit(gtTumor[1],"/")[[1]])
    gtNormalGT <- unlist(strsplit(gtNormal[1],"/")[[1]])
    fixedAlt <- unlist(strsplit(fixed["ALT"],","))
    allelesToKeep <- (order(as.numeric(gtTumorAD)[-1],decreasing = T)[c(1:4)])
    allelesToKeep <- allelesToKeep[order(allelesToKeep)]
    fixedAlt <- fixedAlt[allelesToKeep]
    gtTumorAD <- gtTumorAD[c(1,allelesToKeep+1)]
    gtNormalAD <- gtNormalAD[c(1,allelesToKeep+1)]
    gtTumorAF <- gtTumorAF[c(allelesToKeep+1)]
    gtNormalAF <- gtNormalAF[c(allelesToKeep+1)]
    #gtTumorGT <- gtTumorGT[c(1,allelesToKeep+1)]
    gtTumorGT <- as.character(c(0:4))
    gtTumor[1] <- paste(gtTumorGT,collapse = "/")
    gtNormal[2] <- paste(gtNormalAD,collapse = ",")
    gtTumor[2] <- paste(gtTumorAD,collapse = ",")
    gtNormal[3] <- paste(gtNormalAF,collapse = ",")
    gtTumor[3] <- paste(gtTumorAF,collapse = ",")
    fixed["ALT"] <- paste(fixedAlt,collapse=",")
    fixed["FILTER"] <- paste0(fixed["FILTER"],";more_than_4_alt_alleles")
    vcf@gt[moreThan4[i],normal] <- paste(gtNormal,collapse = ":")
    vcf@gt[moreThan4[i],tumor] <- paste(gtTumor,collapse = ":")
    vcf@fix[moreThan4[i],] <- fixed
  }
  vcfR::write.vcf(vcf,vcfAlissa)
  gzip(fullVcfFile)
}

