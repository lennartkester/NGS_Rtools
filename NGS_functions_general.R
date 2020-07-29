## to do: ##
## in case of misspelled gene in expression check give error message in stead of crash ## 

packages <- c("zoo","readtext","openxlsx","httr","grid","gridExtra","gridBase","pdftools","BiocManager","vcfR","R.utils","colorspace","umap","kknn","viridis")
if (length(setdiff(packages, rownames(installed.packages()))) > 0){
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
bioPackages <- c("BSgenome","BSgenome.Hsapiens.UCSC.hg38","MutationalPatterns")
if (length(setdiff(bioPackages, rownames(installed.packages()))) > 0){
  BiocManager::install(setdiff(bioPackages, rownames(installed.packages())))  
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
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)
library(MutationalPatterns)
library(vcfR)
library(R.utils)
library(colorspace)
library(kknn)
library(viridis)

baseDirWES <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WES/"
baseDirWTS <-  "G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WTS (RNA-Seq)/"
fileList <- list()
fileList$bmlijst <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Isolaties/Lijst BM genereren.xlsx"
fileList$bslijst <- "G://Diagnostisch Lab/Laboratorium/Histologie/BioSource aanmaken/Lijst BS genereren.xlsx"
fileList$NGSdiagnostiek <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Mol-NGS/NGS_Diagnostiek/NGS-Diagnostiek_v3.0.xlsx"
fileList$ither <- "G://Diagnostisch Lab/Laboratorium/Moleculair/Isolaties/iTHER/Ither overzicht.xlsx"
fileList$fusionBlacklist <- paste0(baseDirWTS,"/QualityControl/20200420_STARfusion_blackList_v2.csv")

getFileList <- function(seqRunDir,rootDir,pattern){
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
    perLine <- strsplit(filetext$text,"\n")[[1]]
    targetFiles <- perLine[grep(pattern,perLine)]
    targetFiles <- paste0("http://",sapply(targetFiles, function(x) strsplit(strsplit(x,"http://")[[1]][2],pattern)[[1]][1]),pattern)
    targetFiles <- targetFiles[order(targetFiles)]
    
    return(list("targetFiles"=unique(targetFiles),"perLine"=perLine))
  }  
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
  wesOverview <- as.data.frame(readxl::read_xlsx(paste(baseDirWES,"Overview WES Runs.xlsx",sep="/")))
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

loadSeqFolders <- function(){
  inputChoicesWTS <- list.dirs("G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WTS (RNA-Seq)",recursive = F,full.names = F)
  inputChoicesWES <- list.dirs("G://Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/WES",recursive = F,full.names = F)
  inputChoices <- unique(c(inputChoicesWES,inputChoicesWTS))
  otherDirs <- c("2018","2019","Backup files","QualityControl","Algemene documenten WES diagnostiek")
  inputChoices <- rev(inputChoices[!(inputChoices %in% otherDirs)])
  inputChoices <- rev(inputChoices[order(inputChoices)])
}



