getWESqcData <- function(){
  samples <- list.dirs(paste0(baseDirWES,"QualityControl/multiQCfiles"),recursive = F,full.names = F)
  
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
    }else{
      print(pairedSamples[i])
    }
  }
  rownames(multiQCdatapaired) <- pairedSamples
  multiQCdatapaired[,"Biomaterial ID pair"] <- sapply(rownames(multiQCdatapaired), function(x) paste(strsplit(x,"_")[[1]][c(1,2)],collapse="_"))
  
  #dups <- multiQCdatapaired$`Biomaterial ID`[duplicated(multiQCdatapaired$`Biomaterial ID pair`)]
  multiQCdatapaired <- multiQCdatapaired[!duplicated(multiQCdatapaired$`Biomaterial ID pair`),]
  multiQCdatapaired$sample1 <- sapply(multiQCdatapaired$`Biomaterial ID pair`,function(x) strsplit(x,"_")[[1]][1])
  multiQCdatapaired$sample2 <- sapply(multiQCdatapaired$`Biomaterial ID pair`,function(x) strsplit(x,"_")[[1]][2])
  multiQCdatapaired <- multiQCdatapaired[,grep("OxoG",colnames(multiQCdatapaired),invert = T)]
  colnames(multiQCdatapaired) <- paste("pair_",colnames(multiQCdatapaired),sep="")
  
  singleSamplesBM <- unique(c(multiQCdatapaired$pair_sample1,multiQCdatapaired$pair_sample2))
  missingSamples <- singleSamplesBM[!(singleSamplesBM %in% sapply(singleSamples,function(x) strsplit(x,"_")[[1]][1]))]
  if(length(missingSamples) > 0){
    sampleDirs <- list.dirs(paste0(rootDir,"QualityControl/multiQCfiles/"),full.names = F)
    for(i in 1:length(missingSamples)){
      singleSamples <- c(singleSamples,sampleDirs[grep(paste0("^",missingSamples[i],"_PMCRZ"),sampleDirs)][1])
    }
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
    }else{
      print(singleSamples[i])
    }
  }
  rownames(multiQCdataSingle) <- sapply(singleSamples,function(x) paste(strsplit(x,"_")[[1]][c(1,2)],collapse="_"))
  multiQCdataSingle[,"Biomaterial ID"] <- sapply(multiQCdataSingle$fastqc_Sample, function(x) strsplit(x,"_")[[1]][1])
  
  #dups <- multiQCdataSingle$`Biomaterial ID`[duplicated(multiQCdataSingle$`Biomaterial ID`)]
  multiQCdataSingle <- multiQCdataSingle[!duplicated(multiQCdataSingle$`Biomaterial ID`),]
  
  tempList <- list()
  for ( i in 1:nrow(multiQCdataSingle)){
    tempList[[i]] <- as.vector(c(multiQCdataSingle[i,],multiQCdatapaired[grep(multiQCdataSingle$`Biomaterial ID`[i],multiQCdatapaired$`pair_Biomaterial ID pair`)[1],]))
  }
  
  multiQCdata <- do.call(rbind,tempList)
  multiQCdata <- as.data.frame(apply(multiQCdata,2,unlist),stringsAsFactors = F)
  rownames(multiQCdata) <- multiQCdata$`Biomaterial ID`
  samples <- unique(c(multiQCdatapaired$pair_sample1,multiQCdatapaired$pair_sample2))
  if (!sum(samples %in% multiQCdata$`Biomaterial ID`) == length(samples)){
    missingSample <- samples[!samples %in% multiQCdata$`Biomaterial ID`]
    mostRecentQC <- list.files(paste(rootDir,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
    multiQCdata <- rbind(multiQCdata,qcdataAll[qcdataAll$`Biomaterial ID` == missingSample,colnames(multiQCdata)])
    
  }
  
  BMlist <- as.data.frame(as.matrix(read.xlsx(fileList$bmlijst)),stringsAsFactors = F)
  BMinfo <- BMlist[BMlist$PMABM %in% samples,c("PMCBS","HIX","PMABS","PMABM","Type.BS","PA-nummer","DIN","Datum.isolatie")]
  
  otherSamples <- samples[!(samples %in% BMinfo$PMABM)]
  otherMatrix <- as.data.frame(matrix(ncol=ncol(BMinfo),nrow=length(otherSamples)))
  colnames(otherMatrix) <- colnames(BMinfo)
  otherMatrix$PMABM <- otherSamples
  otherMatrix$PMABS <- otherSamples
  
  BMinfo <- rbind(BMinfo,otherMatrix)
  
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
  otherBS <- multiQCdata$PMABS[!(multiQCdata$PMABS %in% BSlijst$PMABS)]
  if (length(otherBS) > 0){
    otherBSMat <- matrix(ncol=2,nrow=length(otherBS))
    colnames(otherBSMat) <- colnames(BSlijst)
    otherBSMat[,"PMABS"] <- otherBS
    BSlijst <- rbind(BSlijst,otherBSMat)
  }
  rownames(BSlijst) <- BSlijst$PMABS
  
  multiQCdata <- cbind(multiQCdata,BSlijst[multiQCdata$PMABS,"Tumor.%"])
  colnames(multiQCdata)[ncol(multiQCdata)] <- "tumorPerc"
  multiQCdata$tumorPerc <- as.character(multiQCdata$tumorPerc)
  multiQCdata$uniqueReads <- round(as.numeric(multiQCdata$fastqc_Total.Sequences) * as.numeric(multiQCdata$fastqc_total_deduplicated_percentage) / 100000000,0)
  
  date <- gsub("-","_",Sys.Date())
  write.table(multiQCdata,paste0(baseDirWES,"QualityControl/QualityData/",date,"_QCparameters_metaData_merged.csv"),sep="\t")
  write.xlsx(multiQCdata,paste0(baseDirWES,"QualityControl/QualityData/",date,"_QCparameters_metaData_merged.xlsx"))
  
}


makeWTSPlotsAndDataFrame <- function(fusionProbabilityPlots=F){
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
  
  BMlist <- as.matrix(readtext(fileList$bmlijst))
  datum.isolatie <- convertToDate(BMlist[,"Datum.isolatie"])
  dd <- cbind(as.data.frame(BMlist[,c("PMABS","PMABM","RIN")]),datum.isolatie)
  colnames(dd)[2] <- "Biomaterial ID"
  dd <- dd[!is.na(dd$`Biomaterial ID`),]
  idInBoth <- intersect(mergedData$`Biomaterial ID`,dd$`Biomaterial ID`)
  dd <- dd[dd$`Biomaterial ID` %in% idInBoth,]
  mergedData <- mergedData[idInBoth,]
  rownames(dd)<- dd$`Biomaterial ID`
  dd <- dd[rownames(mergedData),]
  
  mergedData <- cbind(mergedData,dd)
  
  mergedData$rf <- 0
  mergedData$rf[mergedData$fusionGenes != "geen_relevante_fusie"] <- 1
  mergedData$picard_RnaSeqMetrics_CODING_BASES <- as.numeric(mergedData$picard_RnaSeqMetrics_CODING_BASES)
  mergedData$picard_insertSize_MEDIAN_INSERT_SIZE <- as.numeric(mergedData$picard_insertSize_MEDIAN_INSERT_SIZE)
  mergedData$fastqc_Total.Sequences <- as.numeric(mergedData$fastqc_Total.Sequences)
  mergedData$general_stats_Picard..general._mqc.generalstats.picard_general.PCT_MRNA_BASES2 <- as.numeric(mergedData$general_stats_Picard..general._mqc.generalstats.picard_general.PCT_MRNA_BASES2)
  mergedData$general_stats_FastQC_mqc.generalstats.fastqc.percent_duplicates2 <- as.numeric(mergedData$general_stats_FastQC_mqc.generalstats.fastqc.percent_duplicates2)
  mergedData$fastqc_total_deduplicated_percentage <- as.numeric(mergedData$fastqc_total_deduplicated_percentage)
  
  
  NGSlist <- read.xlsx(fileList$NGSdiagnostiek,sheet = "RNA Library prep")
  NGSlist$Opmerking[is.na(NGSlist$Opmerking)] <- "geenopmering"
  NGSlist <- NGSlist[NGSlist$Opmerking != "UDI NextFlex",]
  insertSize <- NGSlist[,c("PMABM","Gem.Fragment.lengte.BA","nM.DNA","Volume.(ul)","totaal.in.oprep")]
  colnames(insertSize) <- c("Biomaterial ID","insertSize_BA","concentration_nM","volume_ul","input_ng")
  insertSize$yield <- as.numeric(insertSize$concentration_nM) * as.numeric(insertSize$volume_ul)
  colnames(insertSize) <- c("Biomaterial ID","insertSize_BA","concentration_nM","volume_ul","input_ng","yield(fmol)")
  dups <- rev(duplicated(rev(insertSize$`Biomaterial ID`)))
  insertSize <- insertSize[!dups,]
  idInBoth <- intersect(mergedData$`Biomaterial ID`,insertSize$`Biomaterial ID`)
  insertSize <- insertSize[insertSize$`Biomaterial ID` %in% idInBoth,]
  rownames(insertSize) <- insertSize$`Biomaterial ID`
  mergedData <- mergedData[mergedData$`Biomaterial ID` %in% idInBoth,]
  
  mergedData2 <- cbind(mergedData,insertSize[rownames(mergedData),])
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
  mergedData2$RIN <- as.numeric(as.character(mergedData2$RIN))
  mergedData2$RIN[mergedData2$RIN == 88] <- 8.8
  
  
  
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
        
 
        dev.off()
        
      }
    }
    
  }
  
  destDir <- paste0(baseDirWTS,"QualityControl/QualityData")
  write.table(mergedData2,paste(destDir,"/",date,"_QCparameters_metaData_merged.csv",sep=""),sep="\t")
  write.xlsx(mergedData2,paste(destDir,"/",date,"_QCparameters_metaData_merged.xlsx",sep=""))
  
}


