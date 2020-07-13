
downloadMutect2vcf <- function(folder){
  vcfFiles <- getFileList(folder,rootDir = baseDirWES,pattern = "_WXS.vcf.gz")$targetFiles
  vcfFiles <- vcfFiles[sapply(vcfFiles,function(x) length(strsplit(x,"_")[[1]])) == 4]
  dir.create(paste0(baseDirWES,folder,"/rawData"),showWarnings = F)
  for ( i in 1:length(vcfFiles)){
    expFileName <- sub("http://files.bioinf.prinsesmaximacentrum.nl/WXS/","",vcfFiles[i])
    destFile <- paste(baseDirWES,folder,"/rawData/",expFileName,sep="")
    if ( file.exists(destFile)){
      message(paste(strsplit(expFileName,"_")[[1]][1] ,"already downloaded"))
    }else{
      message(paste("downloading",expFileName))
      GET(vcfFiles[i], authenticate("lkester", "Dm1mYaiS"),write_disk(destFile,overwrite = T))
    }
  }
  
}

processVcf <- function(folder,vcfFile,VAF005=F){
  fullVcfFile <- paste0(baseDirWES,folder,"/rawData/",vcfFile)
  vcfPASS <- sub(".vcf.gz","_PASS.vcf",fullVcfFile)
  vcfVAF005 <- sub(".vcf.gz","_PASS_VAF005.vcf",fullVcfFile)
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  sample_names <- paste(strsplit(vcfFile,"_")[[1]][c(1,2)],collapse = "_")
  tumorSample <- strsplit(sample_names,"_")[[1]][1]
  if(!VAF005){
    if (!(file.exists(vcfPASS))){
      gunzip(fullVcfFile)
      fullVcfFile <- sub(".gz","",fullVcfFile)
      vcf <- vcfR::read.vcfR(fullVcfFile)
      vcfF <- vcf[vcf@fix[,"FILTER"] == "PASS",]
      vcfF2 <- vcfF[grep("|SNV|",vcfF@fix[,"INFO"],fixed = T) ,]
      vcfR::write.vcf(vcfF2,vcfPASS)
      gzip(fullVcfFile)
    }
    vcfs <- read_vcfs_as_granges(vcfPASS, sample_names, ref_genome)
  }else{
    if (!(file.exists(vcfVAF005))){
      gunzip(fullVcfFile)
      fullVcfFile <- sub(".gz","",fullVcfFile)
      vcf <- vcfR::read.vcfR(fullVcfFile)
      vcfF <- vcf[vcf@fix[,"FILTER"] == "PASS",]
      vcfF2 <- vcfF[grep("|SNV|",vcfF@fix[,"INFO"],fixed = T) ,]
      tumorGT <- sapply(vcfF2@gt[,tumorSample],function(x) strsplit(x,":")[[1]][2])
      tumorGT <- sapply(tumorGT,function(x) strsplit(x,","))
      vcfF3 <- vcfF2[lapply(tumorGT,length) ==2,]
      tumorGT <- tumorGT[lapply(tumorGT,length) ==2]
      tumorGT <- apply(do.call(rbind,tumorGT),2,as.numeric)
      tumorVAF <- tumorGT[,2]/(apply(tumorGT,1,sum))
      vcfF4 <- vcfF3[tumorVAF > 0.05,]
      vcfR::write.vcf(vcfF4,vcfVAF005)
      gzip(fullVcfFile)
    }
    vcfs <- read_vcfs_as_granges(vcfVAF005, sample_names, ref_genome)
    
  }
  
  
  
  type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
  mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
  
  return(mut_mat)
}


getMutationalSignature <- function(mut_mat,sample_names,pdf=F,VAF005=F){
  cancer_signatures = read.table("G:/Diagnostisch Lab/Laboratorium/Moleculair/Patientenuitslagen/NGS_Rtools_dev/refFiles/sigProfiler_exome_SBS_signatures.csv", sep = ",", header = TRUE)
  somaticType <- sapply(c(1:nrow(cancer_signatures)), function(x) paste0(substr(cancer_signatures[x,2],1,1),"[",cancer_signatures[x,1],"]",substr(cancer_signatures[x,2],3,3)))
  new_order = match(row.names(mut_mat), somaticType)
  cancer_signatures = cancer_signatures[as.vector(new_order),]
  row.names(cancer_signatures) = somaticType
  cancer_signatures = as.matrix(cancer_signatures[,3:67])
  
  fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
  select <- which(rowSums(fit_res$contribution) > 1)
  totalMuts <- sum(mut_mat)
  
  layout(mat=matrix(nrow=4,ncol=1,data=c(1:4)),heights = c(0.5,5,1.5,1))
  par(mar=c(0,0,0,0))
  plot(1,1,cex=0,axes=F)
  text(1,1,labels = paste("Mutational Signature",sample_names),cex=1.5)
  par(mar=c(1.1,2.1,0,1.1))
  a <- plot_96_profile(mut_mat, condensed = TRUE)
  frame()
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  vp1 <-plotViewport(c(0,0,0,0))
  print(a,vp=vp1)
  par(mar=c(3.1,2.1,2.1,1.1))
  
  set.seed(1111)
  cols <- sample(rainbow_hcl(ncol(cancer_signatures)))
  if(VAF005){
    barTitle <- paste("COSMIC signature contribution -",totalMuts,"PASS filter SNVs with VAF > 5%")
  }else{
    barTitle <- paste("COSMIC signature contribution -",totalMuts,"PASS filter SNVs")
  }
  barplot(as.matrix(fit_res$contribution[select,1]),col=cols[select],horiz = T,axes=F,main=barTitle)
  
  midSegments <- ((cumsum(fit_res$contribution[,1])-c(0,cumsum(fit_res$contribution[,1][-length(fit_res$contribution[,1])])))/2)+c(0,cumsum(fit_res$contribution[,1][-length(fit_res$contribution[,1])]))
  select2 <- which(rowSums(fit_res$contribution) > 0.05*sum(fit_res$contribution))
  axis(1,at=midSegments[select2],labels = names(midSegments)[select2],lwd=0,lwd.ticks=1,las=3,cex.axis=1.5)
  par(mar=c(0,0,1,0))
  plot(1,1,cex=0,axes=F)
  text(1,1,labels = "See https://cancer.sanger.ac.uk/cosmic/signatures/SBS/ for explanation of the signatures",cex=1.5)
  
  if(pdf){
    dev.off()
  }
  
}
