packages <- c("zoo","readtext","openxlsx","httr","grid","gridExtra","gridBase","pdftools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0){
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
library(readtext)
library(httr)
library(openxlsx)
library(zoo)
library(grid)
library(gridExtra)
library(gridBase)
library(pdftools)
library(shiny)

baseDirWES <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WES/"
baseDirWTS <-  "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WTS (RNA-Seq)/"
fileList <- list()
fileList$bmlijst <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Isolaties/Lijst BM genereren.xlsx"
fileList$bslijst <- "G://Diagnostisch Lab/Laboratorium/Histologie/BioSource aanmaken/Lijst BS genereren.xlsx"
fileList$NGSdiagnostiek <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Mol-NGS/NGS_Diagnostiek/NGS-Diagnostiek_v3.0.xlsx"
fileList$ither <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Isolaties/iTHER/Ither overzicht.xlsx"

makeMetaDataWES <- function(seqRunDir=NULL,rootDir=baseDirWES){
  date <- Sys.Date()
  date <- gsub("-","_",date)
  files <- list.files(paste0(rootDir,seqRunDir),
                      full.names = T,
                      recursive = T,
                      pattern = ".docx")
  files <- files[grep("inks",files)]
  files <- files[grep("\\~\\$",files,invert = T)]
  
  if (identical(files,character(0))){
    return(1)
  }else{
    
    filetext <- readtext(files)
    perLine <- strsplit(filetext$text,"http://|.vcf")[[1]]
    vcfFiles <- perLine[grep("WXS.qci",perLine)]
    
    
    for ( i in 1:length(vcfFiles)){
      vcfFiles[i] <- gsub("HYPERLINK","",vcfFiles[i])
      vcfFiles[i] <- paste0("http://files",strsplit(vcfFiles[i],"files")[[1]][2])
      vcfFiles[i] <- paste0(strsplit(vcfFiles[i],"WXS.qci")[[1]][1],"WXS.qci.vcf")
    }
    
    vcfFiles <- unique(vcfFiles)
    
    for ( i in 1:length(vcfFiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",vcfFiles[i])
      destFile <- paste(rootDir,seqRunDir,"/",fileName,sep="")
      GET(vcfFiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
    }
    
    
    perLine <- strsplit(filetext$text,"http://|zip")[[1]]
    multiQCfiles <- perLine[grep("multiqc_data",perLine)]
    
    
    for ( i in 1:length(multiQCfiles)){
      multiQCfiles[i] <- gsub("HYPERLINK","",multiQCfiles[i])
      multiQCfiles[i] <- paste0("http://files",strsplit(multiQCfiles[i],"files")[[1]][2])
      multiQCfiles[i] <- paste0(strsplit(multiQCfiles[i],"multiqc_data")[[1]][1],"multiqc_data.zip")
    }
    
    multiQCfiles <- unique(multiQCfiles)
    multiQCfiles <- multiQCfiles[grep("multi",multiQCfiles)]
    
    downloadedSamples <- list.files(paste0(rootDir,"QualityControl/multiQCfiles/"),pattern = "zip",full.names = T)
    
    for ( i in 1:length(multiQCfiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",multiQCfiles[i])
      qcFileName <- sub(".zip","",paste(strsplit(fileName,split = "_")[[1]][c(1,3,4)],collapse="_"))
      destFile <- paste(rootDir,"QualityControl/multiQCfiles/",fileName,sep="")
      if ( !(destFile %in% downloadedSamples)){
        GET(multiQCfiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
      }
      destDir <- sub("_WXS.multiqc_data.zip","",destFile)
      if ( !dir.exists(destDir)){
        dir.create(destDir,showWarnings = F)  
        unzip(destFile,exdir = destDir)
      }else if(!dir.exists(paste(destDir,"/",qcFileName,sep=""))){
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
    
    multiQCdata <- do.call(rbind,tempList)
    multiQCdata <- as.data.frame(apply(multiQCdata,2,unlist),stringsAsFactors = F)
    rownames(multiQCdata) <- multiQCdata$`Biomaterial ID`
    samples <- c(multiQCdatapaired$pair_sample1,multiQCdatapaired$pair_sample2)
    if (!sum(samples %in% multiQCdata$`Biomaterial ID`) == length(samples)){
      missingSample <- samples[!samples %in% multiQCdata$`Biomaterial ID`]
      mostRecentQC <- list.files(paste(rootDir,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
      mostRecentQC <- mostRecentQC[length(mostRecentQC)]
      qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
      multiQCdata <- rbind(multiQCdata,qcdataAll[qcdataAll$`Biomaterial ID` == missingSample,colnames(multiQCdata)])
      
    }
    
    BMlist <- as.data.frame(as.matrix(read.xlsx(fileList$bmlijst)),stringsAsFactors = F)
    BMinfo <- BMlist[BMlist$PMABM %in% samples,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer","DIN")]
    rownames(BMinfo) <- BMinfo$PMABM
    
    multiQCdata <- cbind(multiQCdata,BMinfo[multiQCdata$`Biomaterial ID`,])
    
    NGSlist <- read.xlsx(fileList$NGSdiagnostiek,sheet = "DNA library prep")
    libData <- NGSlist[,c("PMABM","Gem.Fragment.lengte.BA","nM.DNA","totaal.in.oprep","DIN","totaal.(ng)")]
    colnames(libData) <- c("Biomaterial ID","insertSize_BA","concentration_nM","input_ng","DIN","yield(ng)")
    libData$DIN <- round(as.numeric(libData$DIN),2)
    libData <- libData[!is.na(libData$`Biomaterial ID`),]
    libData <- libData[libData$`Biomaterial ID` %in% samples,]
    dups <- rev(duplicated(rev(libData$`Biomaterial ID`)))
    libData <- libData[!dups,]
    rownames(libData) <- libData$`Biomaterial ID`
    multiQCdata <- cbind(multiQCdata,libData[multiQCdata$`Biomaterial ID`,])
    
    multiQCdata$insertSize_BA <- as.numeric(multiQCdata$insertSize_BA)
    multiQCdata$concentration_nM <- round(as.numeric(multiQCdata$concentration_nM),2)
    multiQCdata$yield <- round(as.numeric(multiQCdata$yield),0)
    
    BSlijst <- read.xlsx(fileList$bslijst)
    BSlijst <- BSlijst[BSlijst$PMABS %in% multiQCdata$PMABS,c("Tumor.%","PMABS")]
    rownames(BSlijst) <- BSlijst$PMABS
    multiQCdata <- cbind(multiQCdata,BSlijst[multiQCdata$PMABS,"Tumor.%"])
    colnames(multiQCdata)[ncol(multiQCdata)] <- "tumorPerc"
    multiQCdata$tumorPerc <- as.character(multiQCdata$tumorPerc)
    multiQCdata$uniqueReads <- round(as.numeric(multiQCdata$fastqc_Total.Sequences) * as.numeric(multiQCdata$fastqc_total_deduplicated_percentage) / 100000000,0)
    
    itherList <- as.data.frame(read.xlsx(fileList$ither),stringsAsFactors=F)
    itherList <- itherList[!is.na(itherList$HIX),]
    itherList$PMABS.tumor[is.na(itherList$PMABS.tumor)] <- "unknown"
    metaData <- as.data.frame(matrix(nrow=nrow(multiQCdatapaired),ncol=17))
    colnames(metaData) <- c("Skion ID","HiX Nr","PMABS tumor","PMABM tumor","PMABS normal","PMABM normal","T-nummer","Vraagstelling","Tumor perc","UniqueReadsTumor(10^6)","UniqueReadsNormal(10^6)","MeanCoverageTumor","MeanCoverageNormal","novelVariants","ContaminationTumor","ContaminationNormal","TumorNormalCheck")
    for ( i in 1:nrow(metaData)){
      metaData[i,c(1,2,3,4,7,9,10,12,14,15,17)] <- multiQCdata[multiQCdata$`Biomaterial ID` == multiQCdatapaired$pair_sample1[i],c("PMCBS","HIX","PMABS","PMABM","PA-nummer","tumorPerc","uniqueReads","picard_HsMetrics_MEAN_TARGET_COVERAGE","pair_gatk_varianteval_novel_sites","verifybamid_FREEMIX","pair_tumor-normal_comparison_Ratio.as.expected.")]
      metaData[i,c(5,6,11,13,16)] <- multiQCdata[multiQCdata$`Biomaterial ID` == multiQCdatapaired$pair_sample2[i],c("PMABS","PMABM","uniqueReads","picard_HsMetrics_MEAN_TARGET_COVERAGE","verifybamid_FREEMIX")]
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
    
    return(metaData)
  }
}






makeFusionExcelFiles <- function(seqRunDir,rootDir = baseDirWTS){
  date <- Sys.Date()
  date <- gsub("-","_",date)
  files <- list.files(paste0(rootDir,seqRunDir),
                      full.names = T,
                      recursive = T,
                      pattern = ".docx")
  files <- files[grep("inks",files)]
  files <- files[grep("\\~\\$",files,invert = T)]
  
  if (identical(files,character(0))){
    return(1)
  }else{
    
    filetext <- readtext(files)
    perLine <- strsplit(filetext$text,"http://|tsv")[[1]]
    fusionFiles <- perLine[grep("star-fusion_predicted.annotated.filtered",perLine)]
    
    for ( i in 1:length(fusionFiles)){
      fusionFiles[i] <- paste0("http://",fusionFiles[i],"tsv")
    }
    fusionFiles <- unique(fusionFiles)
    samples <- sapply(sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",fusionFiles),function(x) strsplit(x,"_")[[1]][1])
    
    dataList <- list()
    multiQCfiles <- perLine[grep("multiqc_data.zip",perLine)]
    multiQCfiles <- paste("http://",unique(sapply(multiQCfiles,function(x) strsplit(x,"multiqc_data.zip")[[1]][1])),"multiqc_data.zip",sep="")
    for ( i in 1:length(multiQCfiles)){
      fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",multiQCfiles[i])
      qcFileName <- sub(".zip","",paste(strsplit(fileName,split = "_")[[1]][c(1,3,4)],collapse="_"))
      destFile <- paste(baseDirWTS,"QualityControl/multiQCfiles/",fileName,sep="")
      GET(multiQCfiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
      destDir <- sub("_RNA-Seq.multiqc_data.zip","",destFile)
      if ( !dir.exists(destDir)){
        dir.create(destDir,showWarnings = F)  
        unzip(destFile,exdir = destDir)
      }else if(!dir.exists(paste(destDir,"/",qcFileName,sep=""))){
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
    libData <- NGSlist[,c("PMABM","Gem.Fragment.lengte.BA","nM.DNA","Volume.(ul)","totaal.in.oprep","RIN")]
    colnames(libData) <- c("Biomaterial ID","insertSize_BA","concentration_nM","volume_ul","input_ng","RIN")
    libData$yield <- as.numeric(libData$concentration_nM) * as.numeric(libData$volume_ul)
    libData$RIN <- round(as.numeric(libData$RIN),2)
    libData <- libData[!is.na(libData$`Biomaterial ID`),]
    libData <- libData[libData$`Biomaterial ID` %in% samples,]
    dups <- rev(duplicated(rev(libData$`Biomaterial ID`)))
    libData <- libData[!dups,]
    
    rownames(libData) <- libData$`Biomaterial ID`
    
    BMlist <- as.data.frame(as.matrix(read.xlsx(fileList$bmlijst)),stringsAsFactors = F)
    datum.isolatie <- convertToDate(BMlist[,"Datum.isolatie"])
    dd <- cbind(as.data.frame(BMlist[,c("PMABS","PMABM","RIN")]),datum.isolatie)
    colnames(dd)[2] <- "Biomaterial ID"
    dd <- dd[!is.na(dd$`Biomaterial ID`),]
    dd <- dd[dd$`Biomaterial ID` %in% samples,]
    rownames(dd) <- dd$`Biomaterial ID`
    
    libData <- cbind(libData,dd[rownames(libData),],as.numeric(multiQCdata[rownames(libData),"Total.Sequences"]) * as.numeric(multiQCdata[rownames(libData),"total_deduplicated_percentage"])/100,as.numeric(multiQCdata[rownames(libData),"Total.Sequences"]))
    colnames(libData)[c((ncol(libData)-1):ncol(libData))] <- c("unique_reads","total_reads")
    
    blackList <- read.csv(paste0(baseDirWTS,"/QualityControl/20200420_STARfusion_blackList_v2.csv"),sep="\t",stringsAsFactors = F)
    
    samplesDone <- c()
    metaList <- list()
    perLine2 <- strsplit(filetext$text,"\n")[[1]]
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
            if (nchar(curFileData[j,k]) > 32766){
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
      tryCatch(write.xlsx(totalMatrix,paste0(baseDirWTS,seqRunDir,"/",sampleName,"_LKR.xlsx"),colNames = F,overwrite=F),error = function(e) NULL ) 
      samplesDone <- c(samplesDone,sampleName)
      #print(paste("Downloaded sample:",sampleName,"and generated excel file"))
    }
    
    #print("Done!")
    
    samples <- BMlist[BMlist$PMABM %in% samplesDone,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer","RIN")]
    
    BSlijst <- read.xlsx(fileList$bslijst)
    BSlijst <- BSlijst[BSlijst$PMABS %in% samples$PMABS,c("Tumor.%","PMABS")]
    rownames(BSlijst) <- BSlijst$PMABS
    itherList <- as.data.frame(read.xlsx(fileList$ither),stringsAsFactors=F)
    itherList <- itherList[!is.na(itherList$HIX),]
    itherList$PMABS.tumor[is.na(itherList$PMABS.tumor)] <- "unknown"
    vraagStelling <- rep("Diagnostiek",nrow(samples))
    for( i in 1:nrow(samples)){
      if(samples$PMABS[i] %in% itherList$PMABS.tumor){
        vraagStelling[i] <- paste("iTHER",itherList$Ither.nummer[itherList$PMABS.tumor == samples$PMABS[i]])
      }
    }
    samples <- as.data.frame(cbind(samples[,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer")],vraagStelling,samples[,c("RIN")],BSlijst[samples$PMABS,"Tumor.%"]),stringsAsFactors=F)
    
    metaData <- do.call(rbind,metaList)
    samples$yield <- round(metaData[samples$PMABM,"yield"],0)
    samples$uniqueReads <- round(metaData[samples$PMABM,"unique_reads"]/1000000,0)
    samples$input <- metaData[samples$PMABM,"input_ng"]
    colnames(samples) <- c("SKION ID","HIX Nr","BioSource ID","Biomaterial ID","Materiaal","Weefsel#","Vraagstelling","RIN","TumorPercentage","yield(fmol)","uniqueReads(10^6)","input(ng)")
    samples$RIN <- round(as.numeric(as.character(samples$RIN)),1)
    otherSamples <- samplesDone[!(samplesDone %in% BMlist$PMABM)]
    otherMatrix <- as.data.frame(matrix(ncol=ncol(samples),nrow=length(otherSamples)))
    colnames(otherMatrix) <- colnames(samples)
    otherMatrix$`Biomaterial ID` <- otherSamples
    otherMatrix$Vraagstelling[grep("PMRBM",otherMatrix$`Biomaterial ID`)] <- "Research"
    otherMatrix$Vraagstelling[grep("PMOBM",otherMatrix$`Biomaterial ID`)] <- "Organoid"
    
    samples <- rbind(samples,otherMatrix)
    
    tryCatch(write.xlsx(samples,paste0(baseDirWTS,seqRunDir,"/metaData_LKR.xlsx"),colNames = T,overwrite=T),error = function(e) print("metadata file already exists") ) 
    #return(paste("Created excelsheets and metadata for",nrow(samples),"samples"))
    return(samples)
  }
}

makeWTSPlotsAndDataFrame <- function(fusionProbabilityPlots=F){
  setwd(paste0(baseDirWTS,"QualityControl/RstudioQualityControl"))
  
  
  ## exract names of the multiqc files from word documents ##
  
  date <- Sys.Date()
  data <- gsub("-","_",date)
  files <- list.files(baseDirWTS,
                      full.names = T,
                      recursive = T,
                      pattern = ".docx")
  
  ## remove incorrect data files and word document related to fusion genes ##
  files <- files[grep("niet correcte|Analyse|\\~\\$",files,invert = T)]
  
  multiQClist <- list()
  for ( i in 1:length(files)){
    filetext <- readtext(files[i])
    perLine <- strsplit(filetext$text,"http://|###")[[1]]
    multiQClist[[i]] <- perLine[grep("RNA-Seq.multiqc_data.zip",perLine)]
  }
  
  multiQCfiles <- unlist(multiQClist)
  for ( i in 1:length(multiQCfiles)){
    multiQCfiles[i] <- gsub("HYPERLINK","",multiQCfiles[i])
    multiQCfiles[i] <- paste0("http://files",strsplit(multiQCfiles[i],"files")[[1]][2])
    multiQCfiles[i] <- paste0(strsplit(multiQCfiles[i],"zip")[[1]][1],"zip")
  }
  
  multiQCfiles <- unique(multiQCfiles)
  multiQCfiles <- multiQCfiles[grep("multi",multiQCfiles)]
  
  downloadedSamples <- list.files(paste0(baseDirWTS,"QualityControl/multiQCfiles/"),pattern = "zip",full.names = T)
  
  for ( i in 1:length(multiQCfiles)){
    fileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/RNA-Seq/","",multiQCfiles[i])
    qcFileName <- sub(".zip","",paste(strsplit(fileName,split = "_")[[1]][c(1,3,4)],collapse="_"))
    destFile <- paste(baseDirWTS,"QualityControl/multiQCfiles/",fileName,sep="")
    if ( !(destFile %in% downloadedSamples)){
      GET(multiQCfiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
    }
    destDir <- sub("_RNA-Seq.multiqc_data.zip","",destFile)
    if ( !dir.exists(destDir)){
      dir.create(destDir,showWarnings = F)  
      unzip(destFile,exdir = destDir)
    }else if(!dir.exists(paste(destDir,"/",qcFileName,sep=""))){
      unzip(destFile,exdir = destDir)
    }
  }
  
  
  dataTypes <- c("fastqc","general_stats","picard_AlignmentSummaryMetrics","picard_baseContent","picard_dups","picard_gcbias","picard_insertSize","picard_RnaSeqMetrics","picard_varientCalling","star")
  dataTypes2 <- paste0("multiqc_",dataTypes,".txt")
  samples <- list.files(paste0(baseDirWTS,"QualityControl/multiQCfiles/"))
  samples <- samples[grep("zip",samples,invert = T)]
  
  
  ## load data from multiQC files ##
  
  dataList <- list()
  for ( i in 1:length(samples)){
    curSampleDir <- list.dirs(paste0(baseDirWTS,"QualityControl/multiQCfiles/",samples[i]),recursive = T)[2]
    if (is.na(curSampleDir)){
      next
    }
    tempList <- list()
    for ( j in 1:length(dataTypes2)){
      tempList[[j]] <- read.csv(paste0(curSampleDir,"/",dataTypes2[j]),sep="\t",stringsAsFactors = F)
      colnames(tempList[[j]]) <- paste(dataTypes[j],colnames(tempList[[j]]),sep="_")
    }
    dataList[[i]] <- as.data.frame(t(unlist(tempList)))
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
  rownames(multiQCdata) <- sapply(samples,function(x) paste(strsplit(x,"_")[[1]][c(1,2)],collapse="_"))
  multiQCdata[,"Biomaterial ID"] <- sapply(multiQCdata$fastqc_Sample, function(x) strsplit(x,"_")[[1]][1])
  
  dups <- multiQCdata$`Biomaterial ID`[duplicated(multiQCdata$`Biomaterial ID`)]
  multiQCdata <- multiQCdata[!(multiQCdata$`Biomaterial ID` %in% dups),]
  
  ## load data from patient results excel sheet ##
  patientResults <- loadRNAseqOverview()
    patientResults <- as.data.frame(patientResults)
  dups <- as.character(patientResults$`Biomaterial ID`[duplicated(patientResults$`Biomaterial ID`)])
  patientResults <- patientResults[!(patientResults$`Biomaterial ID` %in% dups),]
  idInBoth <- intersect(multiQCdata$`Biomaterial ID`,patientResults$`Biomaterial ID`)
  patientResults <- patientResults[patientResults$`Biomaterial ID` %in% idInBoth,]
  rownames(patientResults) <- patientResults$`Biomaterial ID`
  rownames(multiQCdata) <- multiQCdata$`Biomaterial ID`
  
  multiQCdata <- multiQCdata[rownames(patientResults),]
  
  mergedData <- cbind(multiQCdata,patientResults)
  mergedData$fusionGenes <- as.character(mergedData$`Resultaat RNA seq (relevante)`)
  mergedData$fusionGenes[agrep("KIAA1549--BRAF",mergedData$fusionGenes)] <- "KIAA1549--BRAF" 
  mergedData$fusionGenes[agrep("BCR--ABL",mergedData$fusionGenes)] <- "BCR--ABL1" 
  mergedData$fusionGenes[agrep("ETV6--RUNX",mergedData$fusionGenes)] <- "ETV6--RUNX1" 
  mergedData$fusionGenes[agrep("RUNX1--ETV6",mergedData$fusionGenes)] <- "ETV6--RUNX1" 
  mergedData$fusionGenes[agrep("ETV6--NTRK3",mergedData$fusionGenes)] <- "ETV6--NTRK3" 
  mergedData$fusionGenes[agrep("EWSR1--FLI1",mergedData$fusionGenes)] <- "EWSR1--FLI1" 
  mergedData$fusionGenes[agrep("IGH@--CASC11",mergedData$fusionGenes)] <- "IGH@--MYC"
  mergedData$fusionGenes[agrep("IGH@--MYC",mergedData$fusionGenes)] <- "IGH@--MYC"
  mergedData$fusionGenes[agrep("IGH-@-ext--MYC",mergedData$fusionGenes)] <- "IGH@--MYC"
  mergedData$fusionGenes[agrep("KMT2A--AFDN",mergedData$fusionGenes)] <- "KMT2A--AFDN"
  mergedData$fusionGenes[agrep("KMT2A--AFF1",mergedData$fusionGenes)] <- "KMT2A--AFF1"
  mergedData$fusionGenes[agrep("KMT2A--ELL",mergedData$fusionGenes)] <- "KMT2A--ELL"
  mergedData$fusionGenes[agrep("KMT2A--MLLT1",mergedData$fusionGenes)] <- "KMT2A--MLLT1" 
  mergedData$fusionGenes[agrep("NAB2-STAT6",mergedData$fusionGenes)] <- "NAB2-STAT6" 
  mergedData$fusionGenes[agrep("ZNF384--TAF15",mergedData$fusionGenes)] <- "ZNF384--TAF15" 
  mergedData$fusionGenes[agrep("ZMIZ1--ALB1",mergedData$fusionGenes)] <- "ZMIZ1--ABL1"
  mergedData$fusionGenes[agrep("TTYH3--PDGFA",mergedData$fusionGenes)] <- "TTYH3--PDGFA"
  mergedData$fusionGenes[agrep("SET--NUP214",mergedData$fusionGenes)] <- "SET--NUP214"
  mergedData$fusionGenes[agrep("RUNX1--RUNX1T1",mergedData$fusionGenes)] <- "RUNX1--RUNX1T1"
  mergedData$fusionGenes[agrep("PPP1CB--ALK",mergedData$fusionGenes)] <- "PPP1CB--ALK"
  mergedData$fusionGenes[agrep("PML--RARA",mergedData$fusionGenes)] <- "PML--RARA"
  mergedData$fusionGenes[agrep("RARA--PML",mergedData$fusionGenes)] <- "PML--RARA"
  mergedData$fusionGenes[agrep("PAX3--FOXO1",mergedData$fusionGenes)] <- "PAX3--FOXO1"
  mergedData$fusionGenes[agrep("NUP98--NSD1",mergedData$fusionGenes)] <- "NUP98--NSD1"
  mergedData$fusionGenes[agrep("NPM1--ALK",mergedData$fusionGenes)] <- "NPM1--ALK"
  mergedData$fusionGenes[agrep("NAB2--STAT6",mergedData$fusionGenes)] <- "NAB2--STAT6"
  mergedData$fusionGenes[agrep("Geen relevante",mergedData$fusionGenes)] <- "geen_relevante_fusie" 
  mergedData$fusionGenes[grep("Wordt herhaald|ow input|te laag voor|run 169|zie run 171",mergedData$fusionGenes)] <- "geen_relevante_fusie" 
  
  colnames(mergedData) <- gsub("\r","",colnames(mergedData))
  
  BMlist <- as.matrix(readtext(fileList$bslijst))
  datum.isolatie <- convertToDate(BMlist[,"Datum.isolatie"])
  dd <- cbind(as.data.frame(BMlist[,c("PMABS","PMABM","RIN")]),datum.isolatie)
  colnames(dd)[2] <- "Biomaterial ID"
  dd <- dd[!is.na(dd$`Biomaterial ID`),]
  idInBoth <- intersect(mergedData$`Biomaterial ID`,dd$`Biomaterial ID`)
  dd <- dd[dd$`Biomaterial ID` %in% idInBoth,]
  mergedData <- mergedData[idInBoth,]
  rownames(dd)<- dd$`Biomaterial ID`
  dd <- dd[rownames(mergedData),]
  
  mergedData2 <- cbind(mergedData,dd)
  mergedData2$RIN <- as.numeric(as.character(mergedData2$RIN))
  
  mergedData2$RIN[mergedData2$RIN == 88] <- 8.8
  
  mergedData2$rf <- 0
  mergedData2$rf[mergedData2$fusionGenes != "geen_relevante_fusie"] <- 1
  mergedData2$picard_RnaSeqMetrics_CODING_BASES <- as.numeric(mergedData2$picard_RnaSeqMetrics_CODING_BASES)
  mergedData2$picard_insertSize_MEDIAN_INSERT_SIZE <- as.numeric(mergedData2$picard_insertSize_MEDIAN_INSERT_SIZE)
  mergedData2$fastqc_Total.Sequences <- as.numeric(mergedData2$fastqc_Total.Sequences)
  mergedData2$general_stats_Picard..general._mqc.generalstats.picard_general.PCT_MRNA_BASES2 <- as.numeric(mergedData2$general_stats_Picard..general._mqc.generalstats.picard_general.PCT_MRNA_BASES2)
  mergedData2$general_stats_FastQC_mqc.generalstats.fastqc.percent_duplicates2 <- as.numeric(mergedData2$general_stats_FastQC_mqc.generalstats.fastqc.percent_duplicates2)
  mergedData2$fastqc_total_deduplicated_percentage <- as.numeric(mergedData2$fastqc_total_deduplicated_percentage)
  
  
  NGSlist <- read.xlsx(fileList$NGSdiagnostiek,sheet = "RNA Library prep")
  NGSlist$Opmerking[is.na(NGSlist$Opmerking)] <- "geenopmering"
  NGSlist <- NGSlist[NGSlist$Opmerking != "UDI NextFlex",]
  insertSize <- NGSlist[,c("PMABM","Gem.Fragment.lengte.BA","nM.DNA","Volume.(ul)","totaal.in.oprep")]
  colnames(insertSize) <- c("Biomaterial ID","insertSize_BA","concentration_nM","volume_ul","input_ng")
  insertSize$yield <- as.numeric(insertSize$concentration_nM) * as.numeric(insertSize$volume_ul)
  colnames(insertSize) <- c("Biomaterial ID","insertSize_BA","concentration_nM","volume_ul","input_ng","yield(fmol)")
  dups <- rev(duplicated(rev(insertSize$`Biomaterial ID`)))
  insertSize <- insertSize[!dups,]
  idInBoth <- intersect(mergedData2$`Biomaterial ID`,insertSize$`Biomaterial ID`)
  insertSize <- insertSize[insertSize$`Biomaterial ID` %in% idInBoth,]
  rownames(insertSize) <- insertSize$`Biomaterial ID`
  mergedData2 <- mergedData2[mergedData2$`Biomaterial ID` %in% idInBoth,]
  
  mergedData2 <- cbind(mergedData2,insertSize[rownames(mergedData2),])
  mergedData2$insertSize_BA <- as.numeric(mergedData2$insertSize_BA)
  mergedData2$concentration_nM <- as.numeric(mergedData2$concentration_nM)
  
  mergedData2$tumorTypeSimple <- as.character(mergedData2$`Tumor type simple`)
  mergedData2$uniqueReads <- (mergedData2[,"fastqc_Total.Sequences"] * mergedData2[,"fastqc_total_deduplicated_percentage"])/100000000
  colnames(mergedData2)[ncol(mergedData2)] <- "uniqueReads(10^6)"
  

  BSlijst <- read.xlsx(fileList$bslijst)
  BSlijst <- BSlijst[,c("PMABS","HIX","PMCID","Type","Tumor.%","Necrose%")]
  BSlijst <- BSlijst[!is.na(BSlijst$PMABS),]
  BSlijst <- BSlijst[BSlijst$PMABS %in% mergedData2$PMABS,]
  rownames(BSlijst) <- BSlijst$PMABS
  BSlijst$tumorPercFinal <- NA
  BSlijst$tumorPercFinal[grep("%",BSlijst$`Tumor.%`)] <- BSlijst$`Tumor.%`[grep("%",BSlijst$`Tumor.%`)]
  BSlijst$tumorPercFinal <- gsub("%","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub("<","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub(">","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub(" ","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub("\\?","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub(",","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub("~","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub("[a-zA-Z]","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- sapply(BSlijst$tumorPercFinal, function(x) strsplit(x,"-")[[1]][1])
  BSlijst$tumorPercFinal <- gsub("\\(40\\)","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- as.numeric(BSlijst$tumorPercFinal)/100
  BSlijst$tumorPercFinal[grep("%",BSlijst$`Tumor.%`,invert = T)] <- BSlijst$`Tumor.%`[grep("%",BSlijst$`Tumor.%`,invert = T)]
  BSlijst$tumorPercFinal <- gsub("<","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- gsub("[a-zA-Z]","",BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal <- as.numeric(BSlijst$tumorPercFinal)
  BSlijst$tumorPercFinal[is.na(BSlijst$tumorPercFinal)] <- -1
  BSlijst$tumorPercFinal[BSlijst$tumorPercFinal > 1] <- BSlijst$tumorPercFinal[BSlijst$tumorPercFinal > 1]/100
  BSlijst$tumorPercFinal[BSlijst$tumorPercFinal == -1] <- NA
  
  mergedData2 <- merge(mergedData2,BSlijst,by="PMABS",all.x=T)
  
  
  
  ## diagnotic plots to check RNA quality over time ##
  destDir <- paste0(baseDirWTS,"QualityControl/QualityPlots")
  date <- Sys.Date()
  date <- gsub("-","_",date)
  pdf(paste(destDir,"/",date,"_","RNAseqQuality_checksAnalysis.pdf",sep=""),width=8,height=10)
  layout(mat=matrix(ncol=2,nrow=5,data=c(1:10),byrow = T))
  par(mar=c(1,5,3,2))
  
  
  #1 Number of coding bases
  dataToPlot <- mergedData2[,c("datum.isolatie","picard_RnaSeqMetrics_CODING_BASES","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  dataToPlot[,2] <- dataToPlot[,2]/100000000
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab=expression(paste("Bases in coding regions  ","(10"^"9",")")),xlab="",main="Number of bases in coding regions",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab=expression(paste("Bases in coding regions  ","(10"^"9",")")),xlab="",main="Number of bases in coding regions",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #2 percentage coding bases
  dataToPlot <- mergedData2[,c("datum.isolatie","picard_RnaSeqMetrics_PCT_CODING_BASES","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  dataToPlot[,2] <- as.numeric(dataToPlot[,2])
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="Percentage coding bases sequenced",xlab="",main="Number of bases in coding regions (percentage of total)",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="Percentage coding bases sequenced",xlab="",main="Number of bases in coding regions (percentage of total)",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #3 insert size bioinformatics
  dataToPlot <- mergedData2[,c("datum.isolatie","picard_insertSize_MEDIAN_INSERT_SIZE","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  dataToPlot[dataToPlot[,2] > 1000 ,2] <- 1001
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="Median insert size",xlab="",main="Median insert size",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,2] > 1000,1],dataToPlot[dataToPlot[,2] > 1000,2],pch=20,col='red',cex=1.5)
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="Median insert size",xlab="",main="Median insert size",ylim=c(mi,ma))
  points(last6months[last6months[,2] > 1000,1],last6months[last6months[,2] > 1000,2],pch=20,col='red',cex=1.5)
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #4 percentage duplicates
  dataToPlot <- mergedData2[,c("datum.isolatie","fastqc_total_deduplicated_percentage","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  dataToPlot[dataToPlot[,2] > 1000 ,2] <- 1001
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="deduplicated percentage",xlab="",main="deduplicated percentage",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,2] > 1000,1],dataToPlot[dataToPlot[,2] > 1000,2],pch=20,col='red',cex=1.5)
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="deduplicated percentage",xlab="",main="deduplicated percentage",ylim=c(mi,ma))
  points(last6months[last6months[,2] > 1000,1],last6months[last6months[,2] > 1000,2],pch=20,col='red',cex=1.5)
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  
  #5 number of unique reads
  dataToPlot <- mergedData2[,c("datum.isolatie","fastqc_Total.Sequences","rf","fastqc_total_deduplicated_percentage")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  dataToPlot[,2] <- as.numeric(dataToPlot[,2]) * (as.numeric(dataToPlot[,4])) / 100000000
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab=expression(paste("Nr of unique reads  ","(10"^"6",")")),xlab="",main="Unique reads",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab=expression(paste("Nr of unique reads  ","(10"^"6",")")),xlab="",main="Unique reads",ylim=c(mi,ma))
  points(last6months[last6months[,2] > 1000,1],last6months[last6months[,2] > 1000,2],pch=20,col='red',cex=1.5)
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  dev.off()
  
  
  ## quality checks lab ##
  date <- Sys.Date()
  date <- gsub("-","_",date)
  pdf(paste(destDir,"/",date,"_","RNAseqQuality_checksLab.pdf",sep=""),width=8,height=12)
  layout(mat=matrix(ncol=2,nrow=5,data=c(1:10),byrow = T))
  par(mar=c(1,4,3,3))
  
  #1 RIN
  dataToPlot <- mergedData2[,c("datum.isolatie","RIN","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="RIN",xlab="",main="RIN",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="RIN",xlab="",main="RIN",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #2 insertSize
  dataToPlot <- mergedData2[,c("datum.isolatie","insertSize_BA","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="Insertsize",xlab="",main="Insertsize Bioanalyzer",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="Insertsize",xlab="",main="Insertsize Bioanalyzer",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #3 concentration
  dataToPlot <- mergedData2[,c("datum.isolatie","concentration_nM","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="Concentration (nM)",xlab="",main="Concentration",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="Concentration (nM)",xlab="",main="Concentration",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #4 yield
  dataToPlot <- mergedData2[,c("datum.isolatie","yield(fmol)","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="yield (fmol)",xlab="",main="Total library yield",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="yield (fmol)",xlab="",main="Total library yield",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  #5 total reads
  dataToPlot <- mergedData2[,c("datum.isolatie","fastqc_Total.Sequences","rf")]
  dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
  mi <- min(dataToPlot[,2])-(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  ma <- max(dataToPlot[,2])+(0.05*(max(dataToPlot[,2])-min(dataToPlot[,2])))
  if ( mi < 0) {mi <- 0}
  plot(dataToPlot[,1],dataToPlot[,2],pch=20,ylab="Total reads",xlab="",main="Total number of reads",ylim=c(mi,ma))
  points(dataToPlot[dataToPlot[,3] == 1,1],dataToPlot[dataToPlot[,3] == 1,2],pch=20,col='royalblue1')
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)+rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(dataToPlot[order(dataToPlot[,1]),1],rollmean(dataToPlot[order(dataToPlot[,1]),2],50,fill=NA)-rollapply(dataToPlot[order(dataToPlot[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  last6months <- dataToPlot[dataToPlot[,1] > Sys.Date()-(365.25/2),]
  plot(last6months[,1],last6months[,2],pch=20,ylab="Total reads",xlab="",main="Total number of reads",ylim=c(mi,ma))
  points(last6months[last6months[,3] == 1,1],last6months[last6months[,3] == 1,2],pch=20,col='royalblue1')
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA),col='red',lwd=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)+rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  lines(last6months[order(last6months[,1]),1],rollmean(last6months[order(last6months[,1]),2],50,fill=NA)-rollapply(last6months[order(last6months[,1]),2],50,sd,fill=NA),col='red',lwd=1,lty=2)
  legend("topleft",legend=c("No fusion","Fusion"),col=c("black","royalblue1"),pch=20,bty='n',cex=0.5)
  #legend("left",legend=c("moving average","95% CI moving average"),col=c("red","red"),lty=c(1,2),lwd=c(2,1),bty='n')
  
  dev.off()
  
  
  ## fraction of samples with a detected fusion as a function of a quality aspect cut-off in all samples ##
  date <- Sys.Date()
  date <- gsub("-","_",date)
  
  
  hematoTypes <- c("B-ALL","T-ALL","AML","T-LBL","Anaplastic lymphoma","B-cell lymphoma","Hodgkin lymphoma","JMML","Langerhans cell histiocytosis")
  #hematoTypes <- c("B-ALL")
  
  
  mergedData2$tumorTypeSimple[is.na(mergedData2$tumorTypeSimple)] <- "No HiX diagnosis"
  
  if (fusionProbabilityPlots){
    for ( j in c("allPatients","hemato","solideAndNeuro","ALL")){
      if ( j == "allPatients"){
        mergedData3 <- mergedData2
      }else if (j == "hemato"){
        mergedData3 <- mergedData2[ mergedData2$tumorTypeSimple %in% hematoTypes,]
      }else if (j == "solideAndNeuro"){
        mergedData3 <- mergedData2[!(mergedData2$tumorTypeSimple %in% c(hematoTypes,"Not malignant","No HiX diagnosis")),]
      }else if (j == "ALL"){
        mergedData3 <- mergedData2[mergedData2$tumorTypeSimple == "B-ALL",]
      }
      
      for ( i in c("merged")){
        #for ( i in c("single","merged")){ # if single plots are desired use this line in stead of the line above #
        if (i == "merged"){
          pdf(paste(destDir,"/",date,"_","RNAseqQuality_fusionsDetectedSlidingBin_",j,".pdf",sep=""),width=8,height=12)
          layout(mat=matrix(ncol=2,nrow=3,data=c(1:6),byrow = T))
        }
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","RIN")]
        dataToPlot <- dataToPlot[dataToPlot[,3] != 0,]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
        mi <- quantile(dataToPlot[,3],c(0.05,0.95))[1]
        ma <- quantile(dataToPlot[,3],c(0.05,0.95))[2]
        bw <- 2
        fracFusions <- sapply(seq(mi,ma,0.1),function(x) table(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2])[2]/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,0.1),function(x) {
          #sd(dataToPlot[dataToPlot[,2] > x,3])/sqrt(length(dataToPlot[dataToPlot[,2] > x,3]))
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        samplesRetained <- sapply(seq(mi,ma,0.1),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_RIN_",j,".pdf",sep=""))
        }
        par(mar=c(5,5,3,5))
        yma <- max(fracFusionsSD)*1.1
        plot(seq(mi,ma,0.1),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with RIN ~ x (bin width = ",bw,")",sep=""),xlab="RIN",main=j)
        lines(seq(mi,ma,0.1),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,0.1),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,0.1),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,0.1),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","insertSize_BA")]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
        mi <- quantile(dataToPlot[,3],c(0.05,0.95))[1]
        ma <- quantile(dataToPlot[,3],c(0.05,0.95))[2]
        bw <- 50
        fracFusions <- sapply(seq(mi,ma,1),function(x) table(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2])[2]/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,1),function(x) {
          #sd(dataToPlot[dataToPlot[,2] > x,3])/sqrt(length(dataToPlot[dataToPlot[,2] > x,3]))
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        samplesRetained <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_insertSize_BA_",j,".pdf",sep=""))
        }
        yma <- max(fracFusionsSD)*1.1
        plot(seq(mi,ma,1),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with insert size ~ x (bin width = ",bw,")",sep=""),xlab="Insert size Bioanalyzer",main=j)
        lines(seq(mi,ma,1),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,1),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","concentration_nM")]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
        mi <- quantile(dataToPlot[,3],c(0.05,0.9))[1]
        ma <- quantile(dataToPlot[,3],c(0.05,0.9))[2]
        bw <- 20
        fracFusions <- sapply(seq(mi,ma,1),function(x) table(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2])[2]/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,1),function(x) {
          #sd(dataToPlot[dataToPlot[,2] > x,3])/sqrt(length(dataToPlot[dataToPlot[,2] > x,3]))
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        samplesRetained <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_concentraton_nM_",j,".pdf",sep=""))
        }
        yma <- max(fracFusionsSD)*1.1
        plot(seq(mi,ma,1),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with concentration ~ x (bin width = ",bw,")",sep=""),xlab="concentration (nM)",main=j)
        lines(seq(mi,ma,1),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,1),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","yield(fmol)")]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
        mi <- quantile(dataToPlot[,3],c(0.05,0.9))[1]
        ma <- quantile(dataToPlot[,3],c(0.05,0.9))[2]
        bw <- 1000
        fracFusions <- sapply(seq(mi,ma,10),function(x) table(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2])[2]/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,10),function(x) {
          #sd(dataToPlot[dataToPlot[,2] > x,3])/sqrt(length(dataToPlot[dataToPlot[,2] > x,3]))
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        samplesRetained <- sapply(seq(mi,ma,10),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_yield_",j,".pdf",sep=""))
        }
        yma <- max(fracFusionsSD)*1.1
        plot(seq(mi,ma,10),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with yield ~ x (bin width = ",bw,")",sep=""),xlab="yield (fmol)",main=j)
        lines(seq(mi,ma,10),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,10),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,10),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,10),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","picard_dups_READ_PAIRS_EXAMINED","picard_dups_PERCENT_DUPLICATION","tumorTypeSimple")]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot[,c(1:4)]),]
        dataToPlot[,3] <- as.numeric(dataToPlot[,3])*(1-as.numeric(dataToPlot[,4]))/1000000
        mi <- floor(quantile(dataToPlot[,3],c(0.05,0.95)))[1]
        ma <- floor(quantile(dataToPlot[,3],c(0.05,0.95)))[2]
        bw <- 20
        fracFusions <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2] == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,1),function(x) {
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        tumorPerc <- sapply(seq(mi,ma,1),function(x) mean(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),6],na.rm = T))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_NrofUniqueReadPairs_",j,".pdf",sep=""))  
        }
        yma <- max(fracFusionsSD)*1.1
        samplesRetained <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        plot(seq(mi,ma,1),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with unique read pairs ~ x (bin width = ",bw,")",sep=""),xlab=expression(paste("Nr of unique read pairs  ","(10"^"6",")")),main=j)
        lines(seq(mi,ma,1),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,1),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        dataToPlot <- mergedData3[,c("datum.isolatie","rf","fastqc_Total.Sequences","fastqc_total_deduplicated_percentage","tumorTypeSimple")]
        dataToPlot <- dataToPlot[complete.cases(dataToPlot[,c(1:4)]),]
        dataToPlot[,3] <- as.numeric(dataToPlot[,3])*(as.numeric(dataToPlot[,4]))/100000000
        mi <- floor(quantile(dataToPlot[,3],c(0.05,0.95)))[1]
        ma <- floor(quantile(dataToPlot[,3],c(0.05,0.95)))[2]
        bw <- 20
        fracFusions <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2] == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        fracFusionsSD <- sapply(seq(mi,ma,1),function(x) {
          return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        })
        tumorPerc <- sapply(seq(mi,ma,1),function(x) mean(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),6],na.rm = T))
        if ( i == "single"){
          pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_NrofUniqueReads_",j,".pdf",sep=""))  
        }
        yma <- max(fracFusionsSD)*1.1
        samplesRetained <- sapply(seq(mi,ma,1),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        plot(seq(mi,ma,1),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with unique reads ~ x (bin width = ",bw,")",sep=""),xlab=expression(paste("Nr of unique reads  ","(10"^"6",")")),main=j)
        lines(seq(mi,ma,1),fracFusions,pch=20,lwd=2)
        lines(seq(mi,ma,1),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        lines(seq(mi,ma,1),samplesRetained*yma,lwd=2,col='royalblue1')
        axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        if ( i == "single"){
          dev.off()
        }
        
        # if (j %in% c("allPatients","solideAndNeuro")){
        #   dataToPlot <- mergedData3[,c("datum.isolatie","rf","tumorPercFinal")]
        #   dataToPlot <- dataToPlot[complete.cases(dataToPlot),]
        #   mi <- quantile(dataToPlot[,3],c(0.05,0.9))[1]
        #   ma <- quantile(dataToPlot[,3],c(0.05,0.9))[2]
        #   bw <- 0.2
        #   fracFusions <- sapply(seq(mi,ma,0.01),function(x) table(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2])[2]/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2))))
        #   fracFusionsSD <- sapply(seq(mi,ma,0.01),function(x) {
        #     #sd(dataToPlot[dataToPlot[,2] > x,3])/sqrt(length(dataToPlot[dataToPlot[,2] > x,3]))
        #     return(quantile(sapply(c(1:1000), function(y) sum(sample(dataToPlot[dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)),2],replace = T) == 1)/sum(dataToPlot[,3] > (x-(bw/2)) & dataToPlot[,3] < (x+(bw/2)))),probs = c(0.05,0.95)))
        #   })
        #   samplesRetained <- sapply(seq(mi,ma,0.01),function(x) sum(dataToPlot[,3] > x)/nrow(dataToPlot))
        #   if ( i == "single"){
        #     pdf(paste(destDir,"/",date,"_RNAseqQuality_fusionsDetected_yield_",j,".pdf",sep=""))
        #   }
        #   yma <- max(fracFusionsSD)*1.1
        #   plot(seq(mi,ma,0.01),fracFusions,pch=20,cex=0,ylim=c(0,yma),ylab=paste("Fraction of fusion containing samples with tumor fraction ~ x (bin width = ",bw,")",sep=""),xlab="tumor fraction",main=j)
        #   lines(seq(mi,ma,0.01),fracFusions,pch=20,lwd=2)
        #   lines(seq(mi,ma,0.01),fracFusionsSD[1,],pch=20,lwd=1,lty=2)
        #   lines(seq(mi,ma,0.01),fracFusionsSD[2,],pch=20,lwd=1,lty=2)
        #   lines(seq(mi,ma,0.01),samplesRetained*yma,lwd=2,col='royalblue1')
        #   axis(4,at=seq(0,yma,length.out = 6),labels = seq(0,1,0.2),col='royalblue1',col.axis='royalblue1')
        #   mtext("Fraction of samples retained with cutoff x",side = 4,line = 3,col='royalblue1',cex = 0.7)
        #   
        # }
        
        dev.off()
        
      }
    }
    
  }
  
  destDir <- paste0(baseDirWTS,"QualityControl/QualityData")
  write.table(mergedData2,paste(destDir,"/",date,"_QCparameters_metaData_merged.csv",sep=""),sep="\t")
  write.xlsx(mergedData2,paste(destDir,"/",date,"_QCparameters_metaData_merged.xlsx",sep=""))
  
}


loadRNAseqOverview <- function(folder=NULL,samples=NULL,type="biomaterial"){
  rnaSeqOverview <- as.data.frame(read.xlsx(paste(baseDirWTS,"Overview RNAseq Runs.xlsx",sep="/")))
  colnames(rnaSeqOverview) <- rnaSeqOverview[1,]
  rnaSeqOverview <- rnaSeqOverview[c(5:nrow(rnaSeqOverview)),]
  rnaSeqOverview <- rnaSeqOverview[apply(rnaSeqOverview,1,function(x) sum(is.na(x))) != (ncol(rnaSeqOverview)),]
  colnames(rnaSeqOverview)[1] <- "seqRunID"
  for ( i in 1:nrow(rnaSeqOverview)){
    if ( i == 1){
      currentRun <- rnaSeqOverview[i,1]
    }
    if ( is.na(rnaSeqOverview[i,1])){
      rnaSeqOverview[i,1] <- currentRun
    }else if (!is.na(rnaSeqOverview[i,1])){
      currentRun <- rnaSeqOverview[i,1]
    }
  }
  rnaSeqOverview[,1] <- sapply(rnaSeqOverview[,1],function(x) sub("\n","",x))
  rnaSeqOverview[,1] <- sapply(rnaSeqOverview[,1],function(x) sub("\r","",x))
  if (type == "biomaterial"){
    if(!is.null(folder)){
      if(!is.null(samples)){
        rnaSeqOverview <- rnaSeqOverview[rnaSeqOverview[,1] == folder & rnaSeqOverview$`Biomaterial ID` %in% samples,]
      }else{
        rnaSeqOverview <- rnaSeqOverview[rnaSeqOverview[,1] == folder,]
      }
    }else{
      if(!is.null(samples)){
        rnaSeqOverview <- rnaSeqOverview[rnaSeqOverview$`Biomaterial ID` %in% samples,]
      }
    }
  }else if (type == "biosource"){
    if(!is.null(folder)){
      rnaSeqOverview <- rnaSeqOverview[rnaSeqOverview[,1] == folder & rnaSeqOverview$`BioSource ID` %in% samples,]
    }else{
      rnaSeqOverview <- rnaSeqOverview[rnaSeqOverview$`BioSource ID` %in% samples,]
    }
    
  }
  colnames(rnaSeqOverview) <- sapply(colnames(rnaSeqOverview),function(x) sub("\n","",x))
  colnames(rnaSeqOverview) <- sapply(colnames(rnaSeqOverview),function(x) sub("\r","",x))
  
  rnaSeqOverview <- rnaSeqOverview[!is.na(rnaSeqOverview[,"Biomaterial ID"]),]
  rownames(rnaSeqOverview) <- make.unique(rnaSeqOverview[,"Biomaterial ID"])
  rnaSeqOverview <- rnaSeqOverview[,!is.na(colnames(rnaSeqOverview))]
  #rnaSeqOverview <- rnaSeqOverview[samples,]
  colnames(rnaSeqOverview)[16] <- "RIN"
  colnames(rnaSeqOverview)[17] <- "Yield (fmol)"
  colnames(rnaSeqOverview)[19] <- "input RNA (ng)"
  return(rnaSeqOverview)
}


loadWESOverview <- function(folder=NULL,samples=NULL,type="biomaterial"){
  wesOverview <- as.data.frame(read.xlsx(paste(baseDirWES,"Overview WES Runs.xlsx",sep="/")))
  colnames(wesOverview) <- wesOverview[2,]
  wesOverview <- wesOverview[c(5:nrow(wesOverview)),]
  wesOverview <- wesOverview[apply(wesOverview,1,function(x) sum(is.na(x))) != (ncol(wesOverview)),]
  colnames(wesOverview)[1] <- "seqRunID"
  for ( i in 1:nrow(wesOverview)){
    if ( i == 1){
      currentRun <- wesOverview[i,1]
      currentPat <- wesOverview[i,c(2:8)]
    }
    if ( is.na(wesOverview[i,1])){
      wesOverview[i,1] <- currentRun
    }else if (!is.na(wesOverview[i,1])){
      currentRun <- wesOverview[i,1]
    }
    if (sum(is.na(wesOverview[i,c(2:8)])) > 4){
      wesOverview[i,c(2:8)] <- currentPat
    }else if (sum(is.na(wesOverview[i,c(2:8)])) < 4){
      currentPat <- wesOverview[i,c(2:8)]
    }
  }
  wesOverview[,1] <- sapply(wesOverview[,1],function(x) sub("\n","",x))
  wesOverview[,1] <- sapply(wesOverview[,1],function(x) sub("\r","",x))
  if (type == "biomaterial"){
    if(!is.null(folder)){
      if(!is.null(samples)){
        wesOverview <- wesOverview[wesOverview[,1] == folder & wesOverview$`Biomaterial ID` %in% samples,]
      }else{
        wesOverview <- wesOverview[wesOverview[,1] == folder,]
      }
    }else{
      wesOverview <- wesOverview[wesOverview$`Biomaterial ID` %in% samples,]
    }
  }else if (type == "biosource"){
    if(!is.null(folder)){
      wesOverview <- wesOverview[wesOverview[,1] == folder & wesOverview$`BioSource ID` %in% samples,]
    }else{
      wesOverview <- wesOverview[wesOverview$`BioSource ID` %in% samples,]
    }
    
  }
  colnames(wesOverview) <- sapply(colnames(wesOverview),function(x) sub("\n","",x))
  colnames(wesOverview) <- sapply(colnames(wesOverview),function(x) sub("\r","",x))
  
  wesOverview <- wesOverview[!is.na(wesOverview[,"Biomaterial ID"]),]
  wesOverview <- wesOverview[,!is.na(colnames(wesOverview))]
  #wesOverview <- wesOverview[samples %in% wesOverview$`Biomaterial ID`,]
  wesOverview$VAF <- round(as.numeric(wesOverview$VAF),3)
  return(wesOverview)
}


mergeReports <- function(folder=folder, type=type){
  reportDir <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/NGSuitslagen/"
  if ( type == "WTS" | type == "both"){
    if (!dir.exists(paste(baseDirWTS,folder,sep="/"))){
      stop("Specified folder does not exists in WTS (RNA-Seq) folder")
    }
    WTSreports <- list.files(path=paste(baseDirWTS,folder,sep=""),pattern = "WTSreport.pdf",full.names = T)
    WTSreportsIther <- WTSreports[grep("ITHER",WTSreports)]
    WTSreports <- WTSreports[grep("ITHER",WTSreports,invert = T)]
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWTS,folder,"metaData_LKR.xlsx",sep="/")))
    qcdataRun <- qcdataRun[grep("PMRBM|PMOBM",qcdataRun$Biomaterial.ID,invert=T),]
    samples <- unique(sapply(list.files(path=paste(baseDirWTS,folder,sep=""),pattern = "WTSreport.pdf",full.names = F),function(x) strsplit(x = x,"_")[[1]][1]))
    rownames(qcdataRun) <- qcdataRun$Biomaterial.ID
    samplesBiosource <- qcdataRun[samples,"BioSource.ID"]
    wesOverview <- loadWESOverview(samples = samplesBiosource,type="biosource",folder=NULL)
    wesOverview <- wesOverview[!duplicated(wesOverview$`Biomaterial ID`),]
    for ( i in 1:length(samples)){
      sampleBS <- qcdataRun$BioSource.ID[qcdataRun$Biomaterial.ID == samples[i]]
      sampleHIX <- qcdataRun$HIX.Nr[i]
      sampleHIXdir <- paste0(reportDir,sampleHIX)
      ither <- qcdataRun[samples[i],"Vraagstelling"]
      if(!dir.exists(sampleHIXdir)){
        dir.create(sampleHIXdir)
      }
      WESreport <- NULL
      CNVreport <- NULL
      if (sampleBS %in% wesOverview$`BioSource ID`){
        message(paste(samples[i], sampleBS,"has WES data availble"))
        sampleBM <- wesOverview$`Biomaterial ID`[wesOverview$`BioSource ID` == sampleBS]
        sampleBM <- sampleBM[length(sampleBM)]
        WESreport <- paste(baseDirWES,wesOverview$seqRunID[wesOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_WESreport.pdf",sep="")
        WESreport <- WESreport[length(WESreport)]
        if (!file.exists(WESreport)){
          message(paste("WES report for",sampleBS,sampleBM,"does not exist"))
          WESreport <- NULL
        }
        CNVreport <- list.files(path = paste(baseDirWES,wesOverview$seqRunID[wesOverview$`BioSource ID` == sampleBS],"/",sep=""),pattern = "CNVreport.pdf",full.names = T)
        CNVreport <- CNVreport[grep(sampleBM,CNVreport)]
        if (identical(CNVreport, character(0))){
          message(paste("CNV report for",sampleBS,sampleBM,"does not exist"))
          CNVreport <- NULL
        }
      }else{
        message(paste(samples[i],sampleBS,"does not have WES data availble"))
      }
      if(!is.null(CNVreport)){
        CNVtemp <- sub(".pdf","temp.pdf",CNVreport)
        if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek NO"){
          pdftools::pdf_subset(CNVreport,c(1,2,5),CNVtemp)
        }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek HO"){
          pdftools::pdf_subset(CNVreport,c(1,2,4),CNVtemp)
        }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek ST"){
          pdftools::pdf_subset(CNVreport,c(1,2,3),CNVtemp)
        }else{
          pdftools::pdf_subset(CNVreport,c(1,2,3,4,5),CNVtemp)
        }
        if(is.null(WESreport)){
          pdftools::pdf_combine(input = c(WTSreports[i],CNVtemp),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreports[i]))
          pdftools::pdf_combine(input = c(WTSreports[i],CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
        }
        if(!is.null(WESreport)){
          pdftools::pdf_combine(input = c(WTSreports[i],WESreport,CNVtemp),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreports[i]))
          pdftools::pdf_combine(input = c(WTSreports[i],WESreport,CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
        }
        file.remove(CNVtemp)
        ## print WTSreport and CNV report to pdf ##
      }else{
        pdftools::pdf_combine(input = c(WTSreports[i],WESreport),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreports[i]))
        pdftools::pdf_combine(input = c(WTSreports[i],WESreport),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
      }
      if(grepl("ITHER",ither)){
        WESreport <- NULL
        CNVreport <- NULL
        if (sampleBS %in% wesOverview$`BioSource ID`){
          message(paste(samples[i], sampleBS,"has WES data availble"))
          sampleBM <- wesOverview$`Biomaterial ID`[wesOverview$`BioSource ID` == sampleBS]
          sampleBM <- sampleBM[length(sampleBM)]
          WESreport <- paste(baseDirWES,wesOverview$seqRunID[wesOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_",gsub(" ","_",ither),"_WESreport.pdf",sep="")
          WESreport <- WESreport[length(WESreport)]
          if (!file.exists(WESreport)){
            message(paste("WES report for",sampleBS,sampleBM,"does not exist"))
            WESreport <- NULL
          }
          CNVreport <- list.files(path = paste(baseDirWES,wesOverview$seqRunID[wesOverview$`BioSource ID` == sampleBS],"/",sep=""),pattern = "CNVreport.pdf",full.names = T)
          CNVreport <- CNVreport[grep(sampleBM,CNVreport)]
          if (identical(CNVreport, character(0))){
            message(paste("CNV report for",sampleBS,sampleBM,"does not exist"))
            CNVreport <- NULL
          }
        }else{
          message(paste(samples[i],sampleBS,"does not have WES data availble"))
        }
        if(!is.null(CNVreport)){
          CNVtemp <- sub(".pdf","temp.pdf",CNVreport)
          if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek NO"){
            pdftools::pdf_subset(CNVreport,c(1,2,5),CNVtemp)
          }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek HO"){
            pdftools::pdf_subset(CNVreport,c(1,2,4),CNVtemp)
          }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek ST"){
            pdftools::pdf_subset(CNVreport,c(1,2,3),CNVtemp)
          }else{
            pdftools::pdf_subset(CNVreport,c(1,2,3,4,5),CNVtemp)
          }
          if(is.null(WESreport)){
            pdftools::pdf_combine(input = c(WTSreportsIther[grep(samples[i],WTSreportsIther)],CNVtemp),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreportsIther[grep(samples[i],WTSreportsIther)]))
          }
          if(!is.null(WESreport)){
            pdftools::pdf_combine(input = c(WTSreportsIther[grep(samples[i],WTSreportsIther)],WESreport,CNVtemp),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreportsIther[grep(samples[i],WTSreportsIther)]))
          }
          file.remove(CNVtemp)
          ## print WTSreport and CNV report to pdf ##
        }else{
          pdftools::pdf_combine(input = c(WTSreportsIther[grep(samples[i],WTSreportsIther)],WESreport),output = sub("WTSreport.pdf","NGSreport.pdf",WTSreportsIther[grep(samples[i],WTSreportsIther)]))
        }
      }
    }
    reportFiles <- list.files(path=paste(baseDirWTS,folder,sep=""),pattern = "NGSreport.pdf",full.names = T)
    reportFiles <- reportFiles[grep("ITHER",reportFiles,invert = T)]
    reportFiles <- c(paste0(baseDirWTS,folder,"/",folder,"_QCoverview.pdf"),reportFiles)
    pdftools::pdf_combine(input=reportFiles,output=paste0(baseDirWTS,folder,"/",folder,"_runReport.pdf"))
  }
  if ( type == "WES" | type == "both"){
    if (!dir.exists(paste(baseDirWES,folder,sep="/"))){
      stop("Specified folder does not exists in WES folder")
    }
    WESreports <- list.files(path=paste(baseDirWES,folder,sep=""),pattern = "WESreport.pdf",full.names = T)
    WESreportsIther <- WESreports[grep("ITHER",WESreports)]
    WESreports <- WESreports[grep("ITHER",WESreports,invert = T)]
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
    samples <- unique(sapply(list.files(path=paste(baseDirWES,folder,sep=""),pattern = "WESreport.pdf",full.names = F),function(x) strsplit(x = x,"_")[[1]][1]))
    rownames(qcdataRun) <- qcdataRun$PMABM.tumor
    samplesBiosource <- qcdataRun[samples,"PMABS.tumor"]
    wtsOverview <- loadRNAseqOverview(samples = samplesBiosource,type="biosource",folder=NULL)
    wesOverview <- loadWESOverview(samples = samples,folder = folder)
    wesOverview <- wesOverview[!duplicated(wesOverview$`Biomaterial ID`),]
    for ( i in 1:length(samples)){
      curWesReport <- WESreports[i]
      sampleBS <- samplesBiosource[i]
      sampleHIX <- qcdataRun$HiX.Nr[i]
      sampleHIXdir <- paste0(reportDir,sampleHIX)
      ither <- qcdataRun[samples[i],"Vraagstelling"]
      if(!dir.exists(sampleHIXdir)){
        dir.create(sampleHIXdir)
      }
      WTSreport <- NULL
      if (sampleBS %in% wtsOverview$`BioSource ID`){
        message(paste(samples[i], sampleBS,"has RNAseq data availble"))
        sampleBM <- wtsOverview$`Biomaterial ID`[wtsOverview$`BioSource ID` == sampleBS]
        sampleBM <- sampleBM[length(sampleBM)]
        WTSreport <- paste(baseDirWTS,wtsOverview$seqRunID[wtsOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_WTSreport.pdf",sep="")
        WTSreport <- WTSreport[length(WTSreport)]
        if (!file.exists(WTSreport)){
          message(paste("WTS report for",sampleBS,sampleBM,"does not exist"))
          WTSreport <- NULL
        }
      }else{
        message(paste(samples[i],sampleBS,"does not have RNAseq data availble"))
      }
      CNVreport <- list.files(paste(baseDirWES,folder,sep=""),pattern = "CNVreport",full.names = T)
      CNVreport <- CNVreport[grep(samples[i],CNVreport)][1]
      if(!identical(CNVreport,character(0)) & !is.na(CNVreport)){
        CNVtemp <- sub(".pdf","temp.pdf",CNVreport)
        if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek NO"){
          pdftools::pdf_subset(CNVreport,c(1,2,5),CNVtemp)
        }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek HO"){
          pdftools::pdf_subset(CNVreport,c(1,2,4),CNVtemp)
        }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek ST"){
          pdftools::pdf_subset(CNVreport,c(1,2,3),CNVtemp)
        }else{
          pdftools::pdf_subset(CNVreport,c(1,2,3,4,5),CNVtemp)
        }
        if(is.null(WTSreport)){
          pdftools::pdf_combine(input = c(curWesReport,CNVtemp),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
          pdftools::pdf_combine(input = c(curWesReport,CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
        }
        if(!is.null(WTSreport)){
          pdftools::pdf_combine(input = c(WTSreport,curWesReport,CNVtemp),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
          pdftools::pdf_combine(input = c(WTSreport,curWesReport,CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
        }
        file.remove(CNVtemp)
        ## print WTSreport and CNV report to pdf ##
      }else if(!is.null(WTSreport)){
        pdftools::pdf_combine(input = c(WTSreport,curWesReport),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
        pdftools::pdf_combine(input = c(WTSreport,curWesReport),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
      }
      if(grepl("ITHER",ither)){
        curWesReport <- sub("WESreport.pdf",paste0(sub(" ","_",ither),"_WESreport.pdf"),curWesReport)
        WTSreport <- NULL
        if (sampleBS %in% wtsOverview$`BioSource ID`){
          message(paste(samples[i], sampleBS,"has RNAseq data availble"))
          sampleBM <- wtsOverview$`Biomaterial ID`[wtsOverview$`BioSource ID` == sampleBS]
          sampleBM <- sampleBM[length(sampleBM)]
          WTSreport <- paste(baseDirWTS,wtsOverview$seqRunID[wtsOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_",sub(" ","_",ither),"_WTSreport.pdf",sep="")
          WTSreport <- WTSreport[length(WTSreport)]
          if (!file.exists(WTSreport)){
            message(paste("WTS report for",sampleBS,sampleBM,"does not exist"))
            WTSreport <- NULL
          }
        }else{
          message(paste(samples[i],sampleBS,"does not have RNAseq data availble"))
        }
        CNVreport <- list.files(paste(baseDirWES,folder,sep=""),pattern = "CNVreport",full.names = T)
        CNVreport <- CNVreport[grep(samples[i],CNVreport)][1]
        if(!identical(CNVreport,character(0)) & !is.na(CNVreport)){
          CNVtemp <- sub(".pdf","temp.pdf",CNVreport)
          if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek NO"){
            pdftools::pdf_subset(CNVreport,c(1,2,5),CNVtemp)
          }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek HO"){
            pdftools::pdf_subset(CNVreport,c(1,2,4),CNVtemp)
          }else if(wesOverview$Vraagstellling[wesOverview$`BioSource ID`==sampleBS] == "Diagnostiek ST"){
            pdftools::pdf_subset(CNVreport,c(1,2,3),CNVtemp)
          }else{
            pdftools::pdf_subset(CNVreport,c(1,2,3,4,5),CNVtemp)
          }
          if(is.null(WTSreport)){
            pdftools::pdf_combine(input = c(curWesReport,CNVtemp),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
            pdftools::pdf_combine(input = c(curWesReport,CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
          }
          if(!is.null(WTSreport)){
            pdftools::pdf_combine(input = c(WTSreport,curWesReport,CNVtemp),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
            pdftools::pdf_combine(input = c(WTSreport,curWesReport,CNVtemp),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
          }
          file.remove(CNVtemp)
          ## print WTSreport and CNV report to pdf ##
        }else if(!is.null(WTSreport)){
          pdftools::pdf_combine(input = c(WTSreport,curWesReport),output = sub("WESreport.pdf","NGSreport.pdf",curWesReport))
          pdftools::pdf_combine(input = c(WTSreport,curWesReport),output = paste0(sampleHIXdir,"/",sampleHIX,"_",sampleBS,"_",samples[i],".pdf"))
        }
        
      }
    }
    reportFiles <- list.files(path=paste(baseDirWES,folder,sep=""),pattern = "NGSreport.pdf",full.names = T)
    reportFiles <- reportFiles[grep("ITHER",reportFiles,invert = T)]
    reportFiles <- c(paste0(baseDirWES,folder,"/",folder,"_QCoverview.pdf"),reportFiles)
    pdftools::pdf_combine(input=reportFiles,output=paste0(baseDirWES,folder,"/",folder,"_runReport.pdf"))
  }
}


generateReport <- function(folder=folder, type=type){
  if (!(type %in% c("WTS","WES","both"))){
    stop("Type has to be either WTS or WES")
  }
  if ( type == "WTS" | type == "both"){
    if (!dir.exists(paste(baseDirWTS,folder,sep="/"))){
      stop("Specified folder does not exists in WTS (RNA-Seq) folder")
    }
    message("Start generating WTS reports")
    samples <- list.files(paste(baseDirWTS,folder,sep="/"),pattern = ".xlsx")
    samples <- samples[grep("PMABM",samples)]
    samples <- sapply(samples,function(x) strsplit(x,"_")[[1]][1])
    samples <- unique(samples)
    samples <- gsub("~\\$","",samples)
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWTS,folder,"metaData_LKR.xlsx",sep="/")),stringsAsFactors=F)
    mostRecentQC <- list.files(paste(baseDirWTS,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
    rnaSeqOverview <- loadRNAseqOverview(folder=folder,samples=samples)
    for ( sample in samples ){
      printWTSreport(folder=folder,sample=sample,baseDirWTS = baseDirWTS,qcdataRun = qcdataRun, qcdataAll = qcdataAll, rnaSeqOverview = rnaSeqOverview)
      if(grepl("iTHER",qcdataRun[qcdataRun$Biomaterial.ID == sample,"Vraagstelling"])){
        printWTSreport(folder=folder,sample=sample,baseDirWTS = baseDirWTS,qcdataRun = qcdataRun, qcdataAll = qcdataAll, rnaSeqOverview = rnaSeqOverview,ITHER=T)
      }
      message(paste("made",sample,"WTS report"))
    }
    makeWTSoverviewSlide(folder)
  }
  if (type == "WES" | type == "both"){
    if (!dir.exists(paste(baseDirWES,folder,sep="/"))){
      stop("Specified folder does not exists in WES folder")
    }
    message("Start generating WES reports")
    samples <- list.files(paste(baseDirWES,folder,sep="/"),pattern = ".qci")
    samples <- samples[grep("PMABM",samples)]
    samples <- sapply(samples,function(x) strsplit(x,"_")[[1]][1])
    samples <- unique(samples)
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
    mostRecentQC <- list.files(paste(baseDirWES,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
    wesOverview <- loadWESOverview(folder=folder,samples=samples)
    for ( sample in samples ){
      printWESreport(folder=folder,sample=sample,baseDirWES = baseDirWES,qcdataRun = qcdataRun, qcdataAll = qcdataAll, wesOverview = wesOverview)
      if(grepl("iTHER",qcdataRun[qcdataRun$PMABM.tumor == sample,"Vraagstelling"])){
        printWESreport(folder=folder,sample=sample,baseDirWES = baseDirWES,qcdataRun = qcdataRun, qcdataAll = qcdataAll, wesOverview = wesOverview,ITHER=T)
      }  
      message(paste("made",sample,"WES report"))
    }
    makeWESoverviewSlide(folder)
  }
}



qcBarplotWTS <- function(dataRun,dataAll,sample,variable){
  data <- as.numeric(dataAll[,variable])
  data <- data[!is.na(data)]
  magnitude <- floor(max(log10(data)))
  br <- seq(0,1*max(data),length.out=50)
  ma <- (ceiling(max(br)/(10**(magnitude))))*(10**(magnitude))
  if(ma/(2*(10**(magnitude-1))) < 30){
    br <- seq(0,ma,2*(10**(magnitude-2)))
  }else if(ma/(2*(10**(magnitude-1))) > 30 & ma/(2*(10**(magnitude-1))) < 50){
    br <- seq(0,ma,1*(10**(magnitude-1)))
  }else{
    br <- seq(0,ma,2*(10**(magnitude-1)))
  }
  #br <- (10**(magnitude-2))*round(br/(10**(magnitude-2)),0)
  binnedData <- table(cut(data,breaks = br))/length(data)
  a <- barplot(binnedData,col='lightgrey',ylim=c(0,1.1*max(binnedData)),axes=F,names.arg = F,border=F,space = 0)
  rownames(a) <- br[-1]
  axis(1,at=c(0,a[seq(1,length(a),length.out=5)][-1]), labels = c(0,rownames(a)[seq(1,length(a),length.out=5)][-1]))
  axis(2)
  title(xlab=variable,line=2.5,cex.lab=1.2)
  title(ylab="Frequency",line=2.5,cex.lab=1.2)
  box()
  abline(v=(dataRun[dataRun$`Biomaterial.ID` == sample,variable]/max(br)*max(a)),col='red',lwd=2)
  if(variable == "uniqueReads(10^6)"){
    lines(x=c(a["40",1],a["40",1]),y = c(0,0.85*max(binnedData)),lty=2)
    text(x = a["40",1],y = 0.95*max(binnedData),labels = "QC cut-off")
  }
  if(variable == "yield(fmol)"){
    lines(x=c(a["1000",1],a["1000",1]),y = c(0,0.85*max(binnedData)),lty=2)
    text(x = a["1000",1],y = 0.95*max(binnedData),labels = "QC cut-off")
  }
}


qcBarplotWES <- function(dataRun,dataAll,sample,variable){
  if(!(variable %in% c("novelVariants","UniqueReads","MeanCoverage","Contamination"))){
    stop("variable can be either novelVariants, UniqueReads, MeanCoverage or Contamination")
  }
  if (variable == "novelVariants"){
    data <- as.data.frame(cbind(as.numeric(dataAll[,"pair_gatk_varianteval_novel_sites"]),dataAll$Status))
    data <- data[data[,2] == "Tumor",]
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "UniqueReads"){
    data <- as.data.frame(cbind(round(dataAll$fastqc_Total.Sequences * dataAll$fastqc_total_deduplicated_percentage / 100000000,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "MeanCoverage"){
    data <- as.data.frame(cbind(round(dataAll$picard_HsMetrics_MEAN_TARGET_COVERAGE,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "Contamination"){
    data <- as.data.frame(cbind(round(dataAll$verifybamid_FREEMIX,4),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }
  
  data <- data[apply(data,1,function(x) sum(is.na(x))==0),]
  magnitude <- floor(max(log10(data[,1])))
  br <- seq(0,1*max(data[,1]),length.out=50)
  ma <- (ceiling(max(br)/(10**(magnitude))))*(10**(magnitude))
  if(ma/(2*(10**(magnitude-1))) < 20){
    br <- seq(0,ma,5*(10**(magnitude-2)))
  }else if(ma/(2*(10**(magnitude-1))) > 20 & ma/(2*(10**(magnitude-1))) < 50){
    br <- seq(0,ma,1*(10**(magnitude-1)))
  }else{
    br <- seq(0,ma,2*(10**(magnitude-1)))
  }
  #br <- (10**(magnitude-2))*round(br/(10**(magnitude-2)),0)
  binnedDataTumor <- table(cut(data[data[,2] == "Tumor",1],breaks = br))/sum(data[,2] == "Tumor")
  binnedDataNormal <- table(cut(data[data[,2] == "Normal",1],breaks = br))/sum(data[,2] == "Normal")
  binnedData <- rbind(binnedDataTumor,binnedDataNormal)
  binnedData[is.na(binnedData)] <- 0
  a <- barplot(binnedData,col=c('black','lightgrey'),beside=T,ylim=c(0,1.1*max(binnedData)),axes=F,names.arg = rep(F,ncol(binnedData)),border=F,space=c(0,0),legend.text = NULL)  
  b <- t(t(apply(a,2,mean)))
  rownames(b) <- br[-1]
  axis(1,at=c(0,b[seq(1,length(b),length.out=5)][-1]), labels = c(0,rownames(b)[seq(1,length(b),length.out=5)][-1]))
  axis(2)
  title(xlab=variable,line=2.5,cex.lab=1.2)
  title(ylab="Frequency",line=2.5,cex.lab=1.2)
  
  box()
  if (variable == "UniqueReads"){
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"UniqueReadsTumor(10^6)"])/max(br)*max(b)),col='red',lwd=2)  
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"UniqueReadsNormal(10^6)"])/max(br)*max(b)),col='orange',lwd=2)  
    legend("topright",pch=15,col=c("black","lightgrey","red","orange"),legend=c("Tumor (historic)","Normal (historic)","Tumor","Normal"),bty='n')
  }
  if (variable == "MeanCoverage"){
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"MeanCoverageTumor"])/max(br)*max(b)),col='red',lwd=2)  
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"MeanCoverageNormal"])/max(br)*max(b)),col='orange',lwd=2)
    lines(x=c(b["100",1],b["100",1]),y = c(0,0.85*max(binnedData)),lty=2,col="grey")
    text(x = b["100",1],y = 0.95*max(binnedData),labels = "Normal cut-off")
    lines(x=c(b["200",1],b["200",1]),y = c(0,0.85*max(binnedData)),lty=2,col="black")
    text(x = b["200",1],y = 0.95*max(binnedData),labels = "Tumor cut-off")
    legend("topright",pch=15,col=c("black","lightgrey","red","orange"),legend=c("Tumor (historic)","Normal (historic)","Tumor","Normal"),bty='n')
  }
  if (variable == "novelVariants"){
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"novelVariants"])/max(br)*max(b)),col='red',lwd=2)  
  }
  if (variable == "Contamination"){
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"ContaminationTumor"])/max(br)*max(b)),col='red',lwd=2)  
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"ContaminationNormal"])/max(br)*max(b)),col='orange',lwd=2)  
    legend("topright",pch=15,col=c("black","lightgrey","red","orange"),legend=c("Tumor (historic)","Normal (historic)","Tumor","Normal"),bty='n')
    lines(x=c(b["0.02",1],b["0.02",1]),y = c(0,0.85*max(binnedData)),lty=2,col="black")
    text(x = b["0.02",1],y = 0.95*max(binnedData),labels = "QC cut-off")
  }
}


printWTSreport <- function(folder,sample,baseDirWTS,qcdataRun,qcdataAll,rnaSeqOverview,ITHER=F){
  #  png(paste(baseDirWTS,folder,"/",sample,"_WTSreport.png",sep=""), 12, 7, units="in", type="cairo", res=300, bg="white")
  if(ITHER){
    itherNr <- gsub(" ","_",qcdataRun[qcdataRun$Biomaterial.ID == sample,"Vraagstelling"])
    itherNr <- sub("02-","",itherNr)
    pdf(paste(baseDirWTS,folder,"/",sample,"_",itherNr,"_WTSreport.pdf",sep=""), width = 12, height = 7, bg="white")
  }else{
    pdf(paste(baseDirWTS,folder,"/",sample,"_WTSreport.pdf",sep=""), width = 12, height = 7, bg="white")
  }
  layout(mat = matrix(ncol=2,nrow=5,data=c(1,1,2,3,2,4,2,5,6,6),byrow = T),widths = c(1,1),heights = c(0.3,1,1,1,1))
  par(mar=c(0,0,0,0))
  plot(1,cex=0,xlim=c(0,1),ylim=c(0,1),axes=F)
  if(ITHER){
    text(0.5,0.5,paste("HIX",sample,"- WTS Fusion Analysis"),cex=2)
  }else{
    text(0.5,0.5,paste("HIX",qcdataRun[qcdataRun$Biomaterial.ID == sample,"HIX.Nr"],"-",sample,"- WTS Fusion Analysis"),cex=2)
  }
  par(mar=c(5.1,4.1,2.1,2.1))
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  
  if(ITHER){
    tab <- rbind(folder,
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(1,3:6)]),
                 itherNr,
                 t(rnaSeqOverview[sample,c(13)]),
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(8:12)]))
    rownames(tab) <- c("Seq Run ID","SKION ID","Biosource ID","Biomaterial ID","Material type","T-number","Vraagstelling","Diagnosis","RIN","Tumor Cell %","Yield (fmol)","Unique Reads (10^6)","Input amount (ng)")
  }else{
    tab <- rbind(folder,
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(1:7)]),
                 t(rnaSeqOverview[sample,c(13)]),
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(8:12)]))
    rownames(tab) <- c("Seq Run ID","SKION ID","HIX ID","Biosource ID","Biomaterial ID","Material type","T-number","Vraagstelling","Diagnosis","RIN","Tumor Cell %","Yield (fmol)","Unique Reads (10^6)","Input amount (ng)")
    
  }
  tab <- as.table(tab)
  cols <- matrix("black",nrow=nrow(tab),ncol=ncol(tab))
  rownames(cols) <- rownames(tab)
  if(!is.na(tab["Yield (fmol)",])){
    if(as.numeric(tab["Yield (fmol)",]) < 1000){
      cols["Yield (fmol)",] <- "red"
    }
  }else{
    cols["Yield (fmol)",] <- "red"
  }
  if(!is.na(tab["Unique Reads (10^6)",])){
    if(as.numeric(tab["Unique Reads (10^6)",]) < 40){
      cols["Unique Reads (10^6)",] <- "red"
    }
  }else{
    cols["Unique Reads (10^6)",] <- "red"
  }
  if(!is.na(tab["RIN",])){
    if(as.numeric(tab["RIN",]) < 2){
      cols["RIN",] <- "red"
    }else if(as.numeric(tab["RIN",]) < 5){
      cols["RIN",] <- "orange"
    }
  }else{
    cols["RIN",] <- "red"
  }
#  if(!is.na(tab["Tumor Cell %",])){
#    if(as.numeric(tab["Tumor Cell %",]) < 20){
#      cols["Tumor Cell %",] <- "red"
#    }
#  }else{
#    cols["Tumor Cell %",] <- "red"
#  }
  tt <- ttheme_default(core=list(fg_params = list(col = cols)))
  grob <-  tableGrob(tab,cols=NULL,theme = tt)  
  grid.draw(grob)
  popViewport(3)
  par(mar=c(4,4,1,1))
  qcBarplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "RIN")
  qcBarplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "yield(fmol)")
  qcBarplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "uniqueReads(10^6)")
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  tab <- rnaSeqOverview[sample,c(21,22,24)]
  tab[,3] <- sapply(lapply(tab[,3],strwrap, width=50), paste, collapse="\n")
  grob <-  tableGrob(tab,rows=NULL)  
  grid.draw(grob)
  popViewport(3)
  dev.off()
}

printWESreport <- function(folder,sample,baseDirWES,qcdataRun,qcdataAll,wesOverview,ITHER=F){
  #png(paste(baseDirWES,folder,"/",sample,"_WESreport.png",sep=""), 12, 7, units="in", type="cairo", res=300, bg="white")
  if(ITHER){
    itherNr <- gsub(" ","_",qcdataRun[qcdataRun$PMABM.tumor == sample,"Vraagstelling"])
    itherNr <- sub("02-","",itherNr)
    pdf(paste(baseDirWES,folder,"/",sample,"_",itherNr,"_WESreport.pdf",sep=""), width = 12, height = 7, bg="white")
  }else{
    pdf(paste(baseDirWES,folder,"/",sample,"_WESreport.pdf",sep=""), width = 12, height = 7, bg="white")
  }
  layout(mat = matrix(ncol=2,nrow=5,data=c(1,1,2,3,2,4,2,5,6,6),byrow = T),widths = c(1,1),heights = c(0.3,1,1,1,1.5))
  par(mar=c(0,0,0,0))
  plot(1,cex=0,xlim=c(0,1),ylim=c(0,1),axes=F)
  if(ITHER){
    text(0.5,0.5,paste(sample,"- WES Variant Analysis -",wesOverview[which(wesOverview$`Biomaterial ID` == sample)[1],"Panel"]),cex=2)
  }else{
    text(0.5,0.5,paste("HIX",qcdataRun[qcdataRun$PMABM.tumor == sample,"HiX.Nr"],"-",sample,"- WES Variant Analysis -",wesOverview[which(wesOverview$`Biomaterial ID` == sample)[1],"Panel"]),cex=2)
  }
  par(mar=c(5.1,4.1,2.1,2.1))
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  if(ITHER){
    tab <- rbind(folder,
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(1,3:7)]),
                 itherNr,
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(13)][1]),
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(9,12:17)]))
    rownames(tab) <- c("Seq Run ID","SKION ID","Biosource ID Tumor","Biomaterial ID Tumor","Biosource ID Normal","Biomaterial ID Normal","T-number","Vraagstelling","Diagnosis","Tumor Cell %","Mean Coverage Tumor","Mean Coverage Normal","Novel Variants","Contamination Tumor","Contamination Normal","Tumor Normal match")
  }else{
    tab <- rbind(folder,
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(1:8)]),
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(13)][1]),
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(9,12:17)]))
    rownames(tab) <- c("Seq Run ID","SKION ID","HIX ID","Biosource ID Tumor","Biomaterial ID Tumor","Biosource ID Normal","Biomaterial ID Normal","T-number","Vraagstelling","Diagnosis","Tumor Cell %","Mean Coverage Tumor","Mean Coverage Normal","Novel Variants","Contamination Tumor","Contamination Normal","Tumor Normal match")
  }
  tab <- as.table(tab)
  cols <- matrix("black",nrow=nrow(tab),ncol=ncol(tab))
  rownames(cols) <- rownames(tab)
  if(!is.na(tab["Mean Coverage Tumor",])){
    if(as.numeric(tab["Mean Coverage Tumor",]) < 200){
      cols["Mean Coverage Tumor",] <- "red"
    }
  }else{
    cols["Mean Coverage Tumor",] <- "red"
  }
  if(!is.na(tab["Mean Coverage Normal",])){
    if(as.numeric(tab["Mean Coverage Normal",]) < 100){
      cols["Mean Coverage Normal",] <- "red"
    }
  }else{
    cols["Mean Coverage Normal",] <- "red"
  }

  if(!is.na(tab["Contamination Tumor",])){
    if(as.numeric(tab["Contamination Tumor",]) > 2){
      cols["Contamination Tumor",] <- "red"
    }
  }else{
    cols["Contamination Tumor",] <- "red"
  }
  if(!is.na(tab["Contamination Normal",])){
    if(as.numeric(tab["Contamination Normal",]) > 2){
      cols["Contamination Normal",] <- "red"
    }
  }else{
    cols["Contamination Normal",] <- "red"
  }
#  if(!is.na(tab["Tumor Cell %",])){
#    if(as.numeric(tab["Tumor Cell %",]) < 20){
#      cols["Tumor Cell %",] <- "red"
#    }
#  }else{
#    cols["Tumor Cell %",] <- "red"
#  }
  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8,col = cols)),colhead = list(fg_params=list(cex = 0.8)),rowhead = list(fg_params=list(cex = 0.8)))
  grob <-  tableGrob(tab,cols=NULL,theme = mytheme)  
  grid.draw(grob)
  popViewport(3)
  par(mar=c(4,4,1,1))
  qcBarplotWES(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "MeanCoverage")
  qcBarplotWES(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "novelVariants")
  qcBarplotWES(dataRun = qcdataRun,dataAll = qcdataAll,sample = sample,variable = "Contamination")
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  tab <- wesOverview[wesOverview$`Biomaterial ID` == sample,c(16:19,23)]
  tab[,5] <- sapply(lapply(tab[,5],strwrap, width=50), paste, collapse="\n")
  tab[,1] <- sapply(lapply(tab[,1],strwrap, width=50), paste, collapse="\n")
  grob <-  tableGrob(tab,rows=NULL)  
  grid.draw(grob)
  popViewport(3)
  dev.off()
}



plotExpression <- function(gene=NULL,sample=NULL,tumorType=NULL,refData=refCohort,runData=runCohort,pdf=F){
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
  if (!is.null(tumorType)){
    if (pdf){
      pdf(paste0(baseDirWTS,folder,"/",sample,"_",gene,"_",tumorType,".pdf"), width = 12, height = 7, bg="white")
    }
    layout(mat=matrix(ncol=2,nrow=1,data=c(1,2)))
    par(mar=c(5,4,2,0.5))
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='grey',xlab=paste(gene,"- Counts per million"),ylab="Samples - all tumor types",main=paste(gene,"- All tumor types"))
    tumorTypeSamples <- rownames(metaData)[metaData$`Tumor type simple` == tumorType]
    points(dataGeneOrdered[tumorTypeSamples],which(names(dataGeneOrdered) %in% tumorTypeSamples),pch=20)
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=2)
    legend("bottomright",legend = c(tumorType,sample),pch=20,col=c('black','red'),bty='n')
    
    dataGeneOrdered <- dataGeneOrdered[names(dataGeneOrdered) %in% c(tumorTypeSamples,sample)]
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='black',xlab=paste(gene,"- Counts per million"),ylab=paste("Samples -",tumorType),main=paste(gene,"-",tumorType))
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=2)
    legend("bottomright",legend = sample,pch=20,col='red',bty='n')
    if (pdf){
      dev.off()
    }else{
      par(mar=c(5.1,4.1,2.1,2.1))
    }
  }
  if (is.null(tumorType)){
    if (pdf){
      pdf(paste0(baseDirWTS,folder,"/",sample,"_",gene,".pdf"), width = 12, height = 7, bg="white")
    }
    plot(dataGeneOrdered,c(1:length(dataGeneOrdered)),pch=20,col='grey',xlab=paste(gene,"- Counts per million"),ylab="Samples - all tumor types",main=paste(gene,"- All tumor types"))
    points(dataGeneOrdered[sample],which(names(dataGeneOrdered) == sample),col='red',pch=20,cex=2)
    legend("bottomright",legend = c(sample),pch=20,col=c('red'),bty='n')
    if (pdf){
      dev.off()
    }
  }
}



loadRefData <- function(countSet = "20200424_PMCdiag_RNAseq_counts_60357.csv"){
  baseDir <- paste0(baseDirWTS,"QualityControl/expressionData/")
  refFiles <- list.files(baseDir,pattern = "geneExpressionRefData.rds")
  
  refFileDate <- max(sapply(refFiles, function(x) strsplit(x,"_")[[1]][1]))
  refDataDate <- strsplit(countSet,"_")[[1]][1]
  if (as.numeric(refFileDate) > as.numeric(refDataDate)){
    refData <- readRDS(paste0(baseDir,refFiles[which.max(sapply(refFiles, function(x) strsplit(x,"_")[[1]][1]))]))
    refData$counts <- apply(refData$counts,2,function(x) (x/sum(x))*1000000)
    return(refData)
  }else{
    countData <- read.csv(paste0(baseDir,countSet),sep="\t",stringsAsFactors = F)
    samples <- sapply(colnames(countData)[c(3:ncol(countData))],function(x) strsplit(x,"_")[[1]][1])
    dups <- samples[duplicated(samples)]
    samplesDedup <- samples[!(samples %in% dups)]
    
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
    saveRDS(list("counts" = countDataDedup,"metaData"=tumorFusion),paste0(baseDir,date,"_geneExpressionRefData.rds"))
    
    countDataDedupNorm <- apply(countDataDedup,2,function(x) (x/mean(x))*1000000)
    return(list("counts" = countDataDedupNorm,"metaData"=tumorFusion))
  }
}

checkExpressionData <- function(folder){
  expressionFile <- paste0(baseDirWTS,folder,"/expressionData/",folder,"_counts.csv")
  if(!file.exists(expressionFile)){
    #  print("Downloading expression data")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }
}

loadRunExpressionData <- function(folder){
  expressionFile <- paste0(baseDirWTS,folder,"/expressionData/",folder,"_counts.csv")
  if(!file.exists(expressionFile)){
    message("Expression data not available yet, start downloading...")
    counts <- downloadExpressionData(folder)
    write.table(counts,expressionFile,sep="\t")
  }else{
    counts <- read.csv(expressionFile,sep="\t",stringsAsFactors = F)
  }
  return(apply(counts,2,function(x) (x/sum(x))*1000000))
}

downloadExpressionData <- function(folder){
  files <- list.files(paste0(baseDirWTS,folder),
                      full.names = T,
                      recursive = T,
                      pattern = ".docx")
  files <- files[grep("inks",files)]
  files <- files[grep("\\~\\$",files,invert = T)]
  
  if (length(files) == 0){
    stop("No data links file found, are you sure the directory is correct and contains the word document with the links and that the document has \"data links\" in the file name?")
  }
  
  filetext <- readtext(files)
  perLine <- strsplit(filetext$text,"http://|.txt")[[1]]
  expressionFiles <- perLine[grep(".gene_id.exon.counts",perLine,fixed = T)]
  
  
  for ( i in 1:length(expressionFiles)){
    expressionFiles[i] <- gsub("HYPERLINK","",expressionFiles[i])
    expressionFiles[i] <- paste0("http://files",strsplit(expressionFiles[i],"files")[[1]][2],".txt")
  }
  
  expressionFiles <- unique(expressionFiles)
  expressionFiles <- expressionFiles[grep("RNA-Seq.gene_id.exon.counts.txt",expressionFiles,fixed = T)]
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


qcBoxplotWTS <- function(dataRun,dataAll,samples,variable){
  data <- as.numeric(dataAll[,variable])
  boxplot(data,col='lightgrey',ylab=variable,main=variable,cex.lab=1.2,cex.axis=1.2)
  dataRun <- dataRun[!is.na(dataRun[,variable]),]
  samples <- samples[samples %in% rownames(dataRun)]
  points(jitter(rep(1,nrow(dataRun)),5),dataRun[samples,variable],pch=20,col='red',cex=2)
  if (variable == "yield(fmol)"){
    abline(h=1000,lty=2,col='grey')
  }else if(variable == "uniqueReads(10^6)"){
    abline(h=40,lty=2,col='grey')
  }
}

qcBoxplotWES <- function(dataRun,dataAll,samples,variable){
  if (variable == "novelVariants"){
    data <- as.data.frame(cbind(as.numeric(dataAll[,"pair_gatk_varianteval_novel_sites"]),dataAll$Status))
    data <- data[data[,2] == "Tumor",]
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "UniqueReads"){
    data <- as.data.frame(cbind(round(dataAll$fastqc_Total.Sequences * dataAll$fastqc_total_deduplicated_percentage / 100000000,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "MeanCoverage"){
    data <- as.data.frame(cbind(round(dataAll$picard_HsMetrics_MEAN_TARGET_COVERAGE,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "Contamination"){
    data <- as.data.frame(cbind(round(dataAll$verifybamid_FREEMIX,4),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }
  if(variable == "novelVariants"){
    boxplot(list("tumor"=data[data[,2] == "Tumor",1]),col=c('lightgrey'),ylab=variable,main=variable,cex.lab=1.2,cex.axis=1.2,las=2)
  }else{
    boxplot(list("tumor"=data[data[,2] == "Tumor",1],"normal"=data[data[,2] == "Normal",1]),col=c('darkgrey','lightgrey'),ylab=variable,main=variable,cex.lab=1.2,cex.axis=1.2,las=2)
  }
  #dataRun <- dataRun[!is.na(dataRun[,variable]),]
  #samples <- samples[samples %in% rownames(dataRun)]
  if (variable == "novelVariants"){
    points(jitter(rep(1,nrow(dataRun)),5),dataRun[samples,"novelVariants"],pch=20,col='red',cex=2)
  }else if(variable == "UniqueReads"){
    points(jitter(rep(1,nrow(dataRun)),5),dataRun[samples,"UniqueReadsTumor(10^6)"],pch=20,col='red',cex=2)
    points(jitter(rep(2,nrow(dataRun)),5),dataRun[samples,"UniqueReadsNormal(10^6)"],pch=20,col='orange',cex=2)
  }else if(variable == "MeanCoverage"){
    points(jitter(rep(1,nrow(dataRun)),5),dataRun[samples,"MeanCoverageTumor"],pch=20,col='red',cex=2)
    points(jitter(rep(2,nrow(dataRun)),5),dataRun[samples,"MeanCoverageNormal"],pch=20,col='orange',cex=2)
  }else if(variable == "Contamination"){
    points(jitter(rep(1,nrow(dataRun)),5),dataRun[samples,"ContaminationTumor"],pch=20,col='red',cex=2)
    if ( sum(dataRun[samples,"ContaminationTumor"] > max(data[,1],na.rm=T)) > 0){
      points(jitter(rep(1,sum(dataRun[samples,"ContaminationTumor"] > max(data[,1],na.rm=T))),5),rep(max(data[,1],na.rm=T),sum(dataRun[samples,"ContaminationTumor"] > max(data[,1],na.rm=T))),pch=20,col='red',cex=4)
    }
    points(jitter(rep(2,nrow(dataRun)),5),dataRun[samples,"ContaminationNormal"],pch=20,col='orange',cex=2)
    if ( sum(dataRun[samples,"ContaminationNormal"] > max(data[,1],na.rm=T)) > 0){
      points(jitter(rep(2,sum(dataRun[samples,"ContaminationNormal"] > max(data[,1],na.rm=T))),5),rep(max(data[,1],na.rm=T),sum(dataRun[samples,"ContaminationNormal"] > max(data[,1],na.rm=T))),pch=20,col='orange',cex=4)
    }
  }
}

makeWTSoverviewSlide <- function(folder){
  mostRecentQC <- list.files(paste(baseDirWTS,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
  mostRecentQC <- mostRecentQC[length(mostRecentQC)]
  qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
  qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWTS,folder,"metaData_LKR.xlsx",sep="/")))
  qcdataRun <- qcdataRun[grep("PMOBM|PMRBM",qcdataRun$Biomaterial.ID,invert = T),]
  rownames(qcdataRun) <- qcdataRun$Biomaterial.ID
  wtsOverview <- loadRNAseqOverview(folder=folder)
  pdf(paste0(baseDirWTS,folder,"/",folder,"_QCoverview.pdf"),width = 12,height = 7)
  layout(mat = matrix(nrow=2,ncol=4,data=c(1,1,1,1,2,3,4,5),byrow = T),widths = c(3,1,1,1), heights = c(1,5))
  par(mar=c(0,0,0,0))
  plot(1,cex=0,xlim=c(0,1),ylim=c(0,1),axes=F)
  text(0.5,0.5,paste("RNA seq QC overview",folder),cex=2)
  
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  
  tab <- cbind(qcdataRun[,c("HIX.Nr","Vraagstelling")],wtsOverview[rownames(qcdataRun),c("Tumor type")])
  colnames(tab)[3] <- "Tumor type"
  tab[,3] <- sapply(lapply(tab[,3],strwrap, width=40), paste, collapse="\n")
  rownames(tab) <- NULL
  grob <-  tableGrob(tab,rows=NULL)  
  grid.draw(grob)
  popViewport(3)
  par(mar=c(4,4,2,1))  
  qcBoxplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "RIN")
  qcBoxplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "yield(fmol)")
  qcBoxplotWTS(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "uniqueReads(10^6)")
  dev.off()
  
}


makeWESoverviewSlide <- function(folder){
  mostRecentQC <- list.files(paste(baseDirWES,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
  mostRecentQC <- mostRecentQC[length(mostRecentQC)]
  qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
  qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
  qcdataRun <- qcdataRun[grep("PMOBM|PMRBM",qcdataRun$PMABM.tumor,invert = T),]
  rownames(qcdataRun) <- qcdataRun$PMABM.tumor
  wesOverview <- loadWESOverview(folder=folder)
  wesOverview <- wesOverview[!is.na(wesOverview$`Tumor type`),]
  pdf(paste0(baseDirWES,folder,"/",folder,"_QCoverview.pdf"),width = 12,height = 7)
  layout(mat = matrix(nrow=2,ncol=4,data=c(1,1,1,1,2,3,4,5),byrow = T),widths = c(3,1,1,1), heights = c(1,5))
  par(mar=c(0,0,0,0))
  plot(1,cex=0,xlim=c(0,1),ylim=c(0,1),axes=F)
  text(0.5,0.5,paste("WES QC overview",folder),cex=2)
  
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  
  tab <- cbind(qcdataRun[,c("HiX.Nr","Vraagstelling")],wesOverview[rownames(qcdataRun),c("Tumor type")])
  colnames(tab)[3] <- "Tumor type"
  rownames(tab) <- NULL
  grob <-  tableGrob(tab,rows=NULL)  
  grid.draw(grob)
  popViewport(3)
  par(mar=c(4,4,2,1))  
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "MeanCoverage")
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "Contamination")
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "novelVariants")
  dev.off()
  
}

loadSeqFolders <- function(){
  inputChoicesWTS <- list.dirs("G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WTS (RNA-Seq)",recursive = F,full.names = F)
  inputChoicesWES <- list.dirs("G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WES",recursive = F,full.names = F)
  inputChoices <- unique(c(inputChoicesWES,inputChoicesWTS))
  otherDirs <- c("2018","2019","Backup files","QualityControl","Algemene documenten WES diagnostiek")
  inputChoices <- rev(inputChoices[!(inputChoices %in% otherDirs)])
  inputChoices <- rev(inputChoices[order(inputChoices)])
}
