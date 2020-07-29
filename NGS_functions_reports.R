mergeReports <- function(folder=folder, type="both"){
  reportDir <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/NGSuitslagen/"
  if ( type == "WTS" | type == "both"){
    if (!dir.exists(paste(baseDirWTS,folder,sep="/"))){
      stop("Specified folder does not exists in WTS (RNA-Seq) folder")
    }
    WTSreports <- list.files(path=paste(baseDirWTS,folder,sep=""),pattern = "WTSreport.pdf",full.names = T)
    if (identical(WTSreports,character(0))){
      return(paste("No reports present in",folder))
    }
    WTSreportsIther <- WTSreports[grep("iTHER|ITHER",WTSreports)]
    WTSreports <- WTSreports[grep("iTHER|ITHER",WTSreports,invert = T)]
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWTS,folder,"metaData_LKR.xlsx",sep="/")))
    qcdataRun <- qcdataRun[grep("PMOBM",qcdataRun$Biomaterial.ID,invert=T),]
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
      if(grepl("iTHER|ITHER",ither)){
        WESreport <- NULL
        CNVreport <- NULL
        if (sampleBS %in% wesOverview$`BioSource ID`){
          message(paste(samples[i], sampleBS,"has WES data availble"))
          sampleBM <- wesOverview$`Biomaterial ID`[wesOverview$`BioSource ID` == sampleBS]
          sampleBM <- sampleBM[length(sampleBM)]
          WESreport <- paste(baseDirWES,wesOverview$seqRunID[wesOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_",gsub(" 02-","_",ither),"_WESreport.pdf",sep="")
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
    reportFiles <- reportFiles[grep("iTHER|ITHER",reportFiles,invert = T)]
    reportFiles <- c(paste0(baseDirWTS,folder,"/",folder,"_QCoverview.pdf"),reportFiles)
    pdftools::pdf_combine(input=reportFiles,output=paste0(baseDirWTS,folder,"/",folder,"_runReport_WTS.pdf"))
  }
  if ( type == "WES" | type == "both"){
    if (!dir.exists(paste(baseDirWES,folder,sep="/"))){
      stop("Specified folder does not exists in WES folder")
    }
    WESreports <- list.files(path=paste(baseDirWES,folder,sep=""),pattern = "WESreport.pdf",full.names = T)
    WESreportsIther <- WESreports[grep("iTHER|ITHER",WESreports)]
    WESreports <- WESreports[grep("iTHER|ITHER",WESreports,invert = T)]
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
    samples <- unique(sapply(list.files(path=paste(baseDirWES,folder,sep=""),pattern = "WESreport.pdf",full.names = F),function(x) strsplit(x = x,"_")[[1]][1]))
    rownames(qcdataRun) <- qcdataRun$PMABM.tumor
    samplesBiosource <- qcdataRun[samples,"PMABS.tumor"]
    wtsOverview <- loadRNAseqOverview(samples = samplesBiosource,type="biosource",folder=NULL)
    wesOverview <- loadWESOverview(samples = samples,folder = folder)
    wesOverview <- wesOverview[!duplicated(wesOverview$`Biomaterial ID`),]
    WESonlyReports <- c()
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
      if (!is.null(WTSreport) && !grepl(folder,WTSreport)){
        WESonlyReports <- c(WESonlyReports,sub("WESreport.pdf","NGSreport.pdf",curWesReport))
      }
      
      if(grepl("iTHER|ITHER",ither)){
        curWesReport <- sub("WESreport.pdf",paste0(sub(" 02-","_",ither),"_WESreport.pdf"),curWesReport)
        WTSreport <- NULL
        if (sampleBS %in% wtsOverview$`BioSource ID`){
          message(paste(samples[i], sampleBS,"has RNAseq data availble"))
          sampleBM <- wtsOverview$`Biomaterial ID`[wtsOverview$`BioSource ID` == sampleBS]
          sampleBM <- sampleBM[length(sampleBM)]
          WTSreport <- paste(baseDirWTS,wtsOverview$seqRunID[wtsOverview$`BioSource ID` == sampleBS],"/",sampleBM,"_",sub(" 02-","_",ither),"_WTSreport.pdf",sep="")
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
    reportFilesWTS <- list.files(path=paste(baseDirWTS,folder,sep=""),pattern = "NGSreport.pdf",full.names = T)
    reportFilesWTS <- c(reportFilesWTS,WESonlyReports)
    reportFilesWTS <- reportFilesWTS[grep("iTHER|ITHER",reportFilesWTS,invert = T)]
    reportFilesWTS <- c(paste0(baseDirWTS,folder,"/",folder,"_QCoverview.pdf"),paste0(baseDirWES,folder,"/",folder,"_QCoverview.pdf"),reportFilesWTS)
    pdftools::pdf_combine(input=reportFilesWTS,output=paste0("G:/Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/meeting/",folder,"_runReport.pdf"))
    
    reportFiles <- list.files(path=paste(baseDirWES,folder,sep=""),pattern = "NGSreport.pdf",full.names = T)
    reportFiles <- reportFiles[grep("iTHER|ITHER",reportFiles,invert = T)]
    reportFiles <- c(paste0(baseDirWES,folder,"/",folder,"_QCoverview.pdf"),reportFiles)
    pdftools::pdf_combine(input=reportFiles,output=paste0(baseDirWES,folder,"/",folder,"_runReport_WES.pdf"))
  }
  return("Succesfully merged reports")
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
    samples <- samples[grep("PMABM|PMGBM|PMLBM|PMRBM",samples)]
    samples <- sapply(samples,function(x) strsplit(x,"_")[[1]][1])
    samples <- unique(samples)
    samples <- gsub("~\\$","",samples)
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWTS,folder,"metaData_LKR.xlsx",sep="/")),stringsAsFactors=F)
    mostRecentQC <- list.files(paste(baseDirWTS,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
    rnaSeqOverview <- loadRNAseqOverview(folder=folder,samples=samples)
    
    for ( sample in samples ){
      out <- tryCatch(printWTSreport(folder=folder,sample=sample,baseDirWTS = baseDirWTS,qcdataRun = qcdataRun, qcdataAll = qcdataAll, rnaSeqOverview = rnaSeqOverview),error=function(e) return(paste(sample," file not writeable, close file and regenerate reports")))
      if(grepl("iTHER",qcdataRun[qcdataRun$Biomaterial.ID == sample,"Vraagstelling"])){
        out <- tryCatch(printWTSreport(folder=folder,sample=sample,baseDirWTS = baseDirWTS,qcdataRun = qcdataRun, qcdataAll = qcdataAll, rnaSeqOverview = rnaSeqOverview,ITHER=T),error=function(e) return(paste(sample," file not writeable, close file and regenerate reports")))
      }
      if(out == paste(sample," file not writeable")){
        return(out)
      }
      message(paste("made",sample,"WTS report"))
    }
    out <- tryCatch(makeWTSoverviewSlide(folder),error=function(e) return(paste("Overviewslide file not writeable, close file and regenerate reports")))
    return(out)
  }
  if (type == "WES" | type == "both"){
    if (!dir.exists(paste(baseDirWES,folder,sep="/"))){
      stop("Specified folder does not exists in WES folder")
    }
    message("Start generating WES reports")
    samples <- list.files(paste(baseDirWES,folder,sep="/"),pattern = ".qci")
    samples <- samples[grep("PMABM|PMGBM|PMLBM|PMRBM",samples)]
    samples <- sapply(samples,function(x) strsplit(x,"_")[[1]][1])
    samples <- unique(samples)
    qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
    mostRecentQC <- list.files(paste(baseDirWES,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
    mostRecentQC <- mostRecentQC[length(mostRecentQC)]
    qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
    wesOverview <- loadWESOverview(folder=folder,samples=samples)
    for ( sample in samples ){
      out <- tryCatch(printWESreport(folder=folder,sample=sample,baseDirWES = baseDirWES,qcdataRun = qcdataRun, qcdataAll = qcdataAll, wesOverview = wesOverview),error=function(e) return(paste(sample," file not writeable, close file and regenerate reports")))
      if(grepl("iTHER",qcdataRun[qcdataRun$PMABM.tumor == sample,"Vraagstelling"])){
        out <- tryCatch(printWESreport(folder=folder,sample=sample,baseDirWES = baseDirWES,qcdataRun = qcdataRun, qcdataAll = qcdataAll, wesOverview = wesOverview,ITHER=T),error=function(e) return(paste(sample," file not writeable, close file and regenerate reports")))
      }  
      message(paste("made",sample,"WES report"))
    }
    out <- tryCatch(makeWESoverviewSlide(folder),error=function(e) return(paste("Overviewslide file not writeable, close file and regenerate reports")))
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
  if(!is.na(dataRun[dataRun$`Biomaterial.ID` == sample,variable])){
    abline(v=(as.numeric(dataRun[dataRun$`Biomaterial.ID` == sample,variable])/max(br)*max(a)),col='red',lwd=2)
  }
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
    data[data[,1] > 8000,1] <- 8000
  }else if(variable == "UniqueReads"){
    data <- as.data.frame(cbind(round(dataAll$fastqc_Total.Sequences * dataAll$fastqc_total_deduplicated_percentage / 100000000,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "MeanCoverage"){
    data <- as.data.frame(cbind(round(dataAll$picard_HsMetrics_MEAN_TARGET_COVERAGE,0),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
  }else if(variable == "Contamination"){
    data <- as.data.frame(cbind(round(dataAll$verifybamid_FREEMIX,4),dataAll$Status))
    data[,1] <- as.numeric(as.character(data[,1]))
    data[data[,1] > 0.08,1] <- 0.08
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
    dataRun[as.numeric(dataRun[,"novelVariants"]) > 8000,"novelVariants"] <- 8000
    abline(v=(as.numeric(dataRun[dataRun$PMABM.tumor == sample,"novelVariants"])/max(br)*max(b)),col='red',lwd=2)  
  }
  if (variable == "Contamination"){
    dataRun[as.numeric(dataRun[,"ContaminationTumor"]) > 0.08,"ContaminationTumor"]<- 0.08
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
    tab <- rbind(itherNr,
                 t(rnaSeqOverview[sample,c(13)]),folder,
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(1,3:6)]),
                 t(rnaSeqOverview[sample,c(15)]),
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(9:12)]))
    rownames(tab) <- c("Vraagstelling","Diagnosis","Seq Run ID","SKION ID","Biosource ID","Biomaterial ID","Material type","T-number","Tumor Cell %","RIN","Yield (fmol)","Unique Reads (10^6)","Input amount (ng)")
  }else{
    tab <- rbind(t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(2)]),
                 t(rnaSeqOverview[sample,c(13)]),
                 folder,
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(1,3:7)]),
                 t(rnaSeqOverview[sample,c(15)]),
                 t(qcdataRun[qcdataRun$Biomaterial.ID == sample,c(9:12)]))
    rownames(tab) <- c("HIX ID","Diagnosis","Seq Run ID","SKION ID","Biosource ID","Biomaterial ID","Material type","T-number","Vraagstelling","Tumor Cell %","RIN","Yield (fmol)","Unique Reads (10^6)","Input amount (ng)")
    
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
  tt <- ttheme_default(core=list(fg_params = list(col = cols,fontface=c(rep("bold", 2), rep("plain",nrow(tab)-2)))))
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
  return("Succesfully generated reports")
}

printWESreport <- function(folder,sample,baseDirWES,qcdataRun,qcdataAll,wesOverview,ITHER=F){
  
  if (sum(wesOverview$`Biomaterial ID` == sample) > 4){
    additionalHeight <- (sum(wesOverview$`Biomaterial ID` == sample)-4)/2.5
  }else{
    additionalHeight <- 0
  }
  #png(paste(baseDirWES,folder,"/",sample,"_WESreport.png",sep=""), 12, 7, units="in", type="cairo", res=300, bg="white")
  if(ITHER){
    itherNr <- gsub(" ","_",qcdataRun[qcdataRun$PMABM.tumor == sample,"Vraagstelling"])
    itherNr <- sub("02-","",itherNr)
    pdf(paste(baseDirWES,folder,"/",sample,"_",itherNr,"_WESreport.pdf",sep=""), width = 12, height = (7 + additionalHeight), bg="white")
  }else{
    pdf(paste(baseDirWES,folder,"/",sample,"_WESreport.pdf",sep=""), width = 12, height = (7 + additionalHeight), bg="white")
  }
  layout(mat = matrix(ncol=2,nrow=5,data=c(1,1,2,3,2,4,2,5,6,6),byrow = T),widths = c(1,1),heights = c(c(0.3,1.2,1.2,1.2),(1.5 + additionalHeight)/((7+additionalHeight)/7)))
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
    tab <- rbind(itherNr,
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(13)][1]),
                 folder,
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(1,3:7)]),
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(14)][1]),
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(12:17)]))
    rownames(tab) <- c("Vraagstelling","Diagnosis","Seq Run ID","SKION ID","Biosource ID Tumor","Biomaterial ID Tumor","Biosource ID Normal","Biomaterial ID Normal","T-number","Tumor Cell %","Mean Coverage Tumor","Mean Coverage Normal","Novel Variants","Contamination Tumor","Contamination Normal","Tumor Normal match")
  }else{
    tab <- rbind(t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(2)]),
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(13)][1]),
                 folder,
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(1,3:8)]),
                 t(wesOverview[wesOverview$`Biomaterial ID` == sample,c(14)][1]),
                 t(qcdataRun[qcdataRun$PMABM.tumor == sample,c(12:17)]))
    rownames(tab) <- c("HIX ID","Diagnosis","Seq Run ID","SKION ID","Biosource ID Tumor","Biomaterial ID Tumor","Biosource ID Normal","Biomaterial ID Normal","T-number","Vraagstelling","Tumor Cell %","Mean Coverage Tumor","Mean Coverage Normal","Novel Variants","Contamination Tumor","Contamination Normal","Tumor Normal match")
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
  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8,col = cols,fontface=c(rep("bold", 2), rep("plain",nrow(tab)-2)))),colhead = list(fg_params=list(cex = 0.8)),rowhead = list(fg_params=list(cex = 0.8)))
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
  return("Succesfully generated reports")
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
    boxplot(log10(data[data[,2] == "Tumor",1]),col=c('lightgrey'),ylab=paste("log10",variable),main=variable,cex.lab=1.2,cex.axis=1.2,las=2)
  }else{
    boxplot(list("tumor"=data[data[,2] == "Tumor",1],"normal"=data[data[,2] == "Normal",1]),col=c('darkgrey','lightgrey'),ylab=variable,main=variable,cex.lab=1.2,cex.axis=1.2,las=2)
  }
  #dataRun <- dataRun[!is.na(dataRun[,variable]),]
  #samples <- samples[samples %in% rownames(dataRun)]
  if (variable == "novelVariants"){
    points(jitter(rep(1,nrow(dataRun)),5),log10(as.numeric(dataRun[samples,"novelVariants"])),pch=20,col='red',cex=2)
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
  qcdataRun <- qcdataRun[grep("PMOBM",qcdataRun$Biomaterial.ID,invert = T),]
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
  return("Succesfully generated reports")
}


makeWESoverviewSlide <- function(folder){
  mostRecentQC <- list.files(paste(baseDirWES,"QualityControl/QualityData",sep="/"),pattern = ".csv",full.names = T)
  mostRecentQC <- mostRecentQC[length(mostRecentQC)]
  qcdataAll <- read.csv(mostRecentQC,stringsAsFactors = F,sep="\t",check.names = F)
  qcdataRun <- as.data.frame(read.xlsx(paste(baseDirWES,folder,"metaData_LKR.xlsx",sep="/")))
  qcdataRun <- qcdataRun[grep("PMOBM",qcdataRun$PMABM.tumor,invert = T),]
  rownames(qcdataRun) <- qcdataRun$PMABM.tumor
  wesOverview <- loadWESOverview(folder=folder)
  wesOverview <- wesOverview[!is.na(wesOverview$`Tumor type`),]
  rownames(wesOverview) <- wesOverview$`Biomaterial ID`
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
  tab[,3] <- sapply(lapply(tab[,3],strwrap, width=40), paste, collapse="\n")
  rownames(tab) <- NULL
  grob <-  tableGrob(tab,rows=NULL)  
  grid.draw(grob)
  popViewport(3)
  par(mar=c(4,4,2,1))  
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "MeanCoverage")
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "Contamination")
  qcBoxplotWES(dataRun = qcdataRun,dataAll = qcdataAll,samples = rownames(qcdataRun),variable = "novelVariants")
  dev.off()
  return("Succesfully generated reports")
}
