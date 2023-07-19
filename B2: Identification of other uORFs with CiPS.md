Identification of other uORFs with CiPS

```
########################################
########################################
# Section 2: Find translated >= 2aa uORFs
# Map Ribo-seq reads to each length of uORFs
# From 2AA to more
rm(list=ls())
# library(dplyr)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(Biostrings)
# library(Rsamtools)
# library(ORFik)
# library(RiboPlotR)
# library(future.apply)
# library(ggplot2)
# library(clipr)
# library(openxlsx)
GTF_path <- "~/Desktop/CTRL_TPC/data/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf"
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF(file = GTF_path,format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByTx <- exonsBy(txdb, by="tx", use.names=T)
length(exonByTx) #21845
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTR)
fiveUTR_ORFs <- findMapORFs(fiveUTR, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F)
fiveUTR_gene_id <- substr(names(fiveUTR),1,9)
fiveUTR_ORFs_gene_id <- substr(names(fiveUTR_ORFs),1,9)
fiveUTR_no_ORFs_gene_id <-  setdiff(fiveUTR_gene_id,fiveUTR_ORFs_gene_id)

length(fiveUTR) #20938
length(fiveUTR_ORFs) #28228
uORF_length <- sum(width(fiveUTR_ORFs))
fiveUTR_ORFs_seq <- extractTranscriptSeqs(FA,fiveUTR_ORFs)

Ribo1=read.delim(file="~/Desktop/CTRL_v1/CTRL_expressed_P_sites_sort_count",header=F,stringsAsFactors=F,sep="\t")
colnames(Ribo1) <- c("count","chr","start","strand")

Ribo1$end <- Ribo1$start
head(Ribo1)
RiboR <- makeGRangesFromDataFrame(Ribo1,keep.extra.columns=T)
RiboR

####
XAA <-function(x,y){
  # print(x)
  STRAND=as.character(strand(fiveUTR_ORFs_xaa[[x]])[1])
  if(STRAND=="+") {z=fiveUTR_ORFs_xaa[[x]]} 
  else if (STRAND=="-"){z=rev(fiveUTR_ORFs_xaa[[x]])}
  aa=unlist(Map(`:`, start(z), end(z)))
  gr=GRanges(seqnames=suppressWarnings(as.character(seqnames(z))[1]), ranges=IRanges(start = aa,width=1),strand=STRAND)
  if(STRAND=="+"){
    gr$frame=rep(1:3,time=length(aa)/3)
    gr1 <- gr[1:(3*(y-1))]
    gr2 <- gr[(3*(y-1)+1):(3*y)]
    gr3 <- gr[(3*y+1):(3*(y+1))]
  } else if (STRAND=="-"){
    gr$frame=rep(3:1,time=length(aa)/3)
    gr1 <- gr[(3*(y+1)):7]
    gr2 <- gr[6:4]
    gr3 <- gr[3:1]
  }
  df0 <- data.frame(count=c(0,0,0),frame=c(1,2,3))
  ####***
  ranges1 <- subsetByOverlaps(RiboR,gr1)
  hits1 <- findOverlaps(RiboR, gr1)
  xframe1 <- unlist(IntegerList(split(gr1$frame[subjectHits(hits1)], queryHits(hits1))))
  mcols(ranges1) <- DataFrame(mcols(ranges1), xframe1)
  
  df1 <- data.frame(count=ranges1$count,frame=ranges1$xframe1)
  df1 <- rbind(df1,df0)
  df1t1 <- table(df1$frame) %>% as.data.frame()
  df1t2 <- df1 %>% group_by(frame) %>% 
    summarise(total_counts = sum(count)) %>% as.data.frame()
  df1t2$frame_freq <- df1t1$Freq-1
  
  ####***
  ranges2 <- subsetByOverlaps(RiboR,gr2)
  hits2 <- findOverlaps(RiboR, gr2)
  xframe2 <- unlist(IntegerList(split(gr2$frame[subjectHits(hits2)], queryHits(hits2))))
  mcols(ranges2) <- DataFrame(mcols(ranges2), xframe2)
  df2 <- data.frame(count=ranges2$count,frame=ranges2$xframe2)
  df2 <- rbind(df2,df0)
  df2t1 <- table(df2$frame) %>% as.data.frame()
  df2t2 <- df2 %>% group_by(frame) %>% 
    summarise(total_counts = sum(count)) %>% as.data.frame()
  df2t2$frame_freq <- df2t1$Freq-1
  
  ####***
  ranges3 <- subsetByOverlaps(RiboR,gr3)
  hits3 <- findOverlaps(RiboR, gr3)
  xframe3 <- unlist(IntegerList(split(gr3$frame[subjectHits(hits3)], queryHits(hits3))))
  mcols(ranges3) <- DataFrame(mcols(ranges3), xframe3)
  df3 <- data.frame(count=ranges3$count,frame=ranges3$xframe3)
  df3 <- rbind(df3,df0)
  df3t1 <- table(df3$frame) %>% as.data.frame()
  df3t2 <- df3 %>% group_by(frame) %>% 
    summarise(total_counts = sum(count)) %>% as.data.frame()
  df3t2$frame_freq <- df3t1$Freq-1
  
  DF <- cbind(df1t2,df2t2[,2:3],df3t2[,2:3])
  colnames(DF) <- c("Frame","Counts1","FrameFreq1","CountsLast","FrameFreqLast","CountStop","FrameFreqStop")
  DF
}

#dir.create("~/Desktop/CTRL_TPC/data/Rdata")
Length <- as.data.frame(table(uORF_length),stringsAsFactors = FALSE)
Length$aa <- (as.numeric(Length$uORF_length)/3)-1
length(Length$aa) #139
AA <- Length$aa[-1][1:(length(Length$aa)-1)] #-1 to remove minimum uORFs
plan(multisession)
for(i in AA){
  gc()
  print(paste0("AA= ",i))
  fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length==(i+1)*3]
  print(length(fiveUTR_ORFs_xaa))
  #future_lapply
  B=future_lapply(seq_len(length(fiveUTR_ORFs_xaa)),function(x) XAA(x,y=i))
  names(B) <- names(fiveUTR_ORFs_xaa)
  save(B,file=paste0("~/Desktop/CTRL_TPC/data/Rdata/Mar02_2023/2Max_isoform_0.5TPM_AA_",i,"_.RData"))
}

####################################################################
### Functions to get counts from different frames or total counts
InFrameCounts <- function(x,A){
  N=A[x]
  ToTal_in_frame_counts_wo_stop <- sum(N[[1]][1,2],N[[1]][1,4],N[[1]][2,4])
}

TotalCounts <- function(x,A){
  N=A[x]
  ToTal_counts_wo_stop <-sum(N[[1]][,c(2,4)])
}

### Functions to get in frames sites with Ribo-seq reads or total sites
InFrameSites <- function(x,A){
  N=A[x]
  ToTal_in_frame_sites_wo_stop <- sum(N[[1]][1,3],N[[1]][1,5],N[[1]][2,5])
}

TotalSites <- function(x,A){
  N=A[x]
  ToTal_sites_wo_stop <-sum(N[[1]][,c(3,5)])
}

InFrameCountSites <- function(z){
  print(z)
  fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length==(z+1)*3]
  load(paste0("~/Desktop/CTRL_TPC/data/Rdata/2Max_isoform_0.5TPM_AA_",z,"_.RData")) #load B
  InFrameC_AA <<- unlist(lapply(seq_len(length(B)), function(x) InFrameCounts(x,A=B)))
  TotalC_AA <<- unlist(lapply(seq_len(length(B)), function(x) TotalCounts(x,A=B)))
  
  InFrameS_AA <<- unlist(lapply(seq_len(length(B)), function(x) InFrameSites(x,A=B)))
  TotalS_AA <<- unlist(lapply(seq_len(length(B)), function(x) TotalSites(x,A=B)))
  STRAND <<- sapply(seq_len(length(B)), function(x) as.character(strand(fiveUTR_ORFs_xaa[names(B)[x]])[[1]][1]))
  INFO <- data.frame(ORF_id =names(B),
                     InFrameC=InFrameC_AA,
                     TotalC=TotalC_AA,
                     InFrameC_perc = round(InFrameC_AA/TotalC_AA*100,1),
                     InFrameS=InFrameS_AA,
                     TotalS=TotalS_AA,
                     InFrameS_perc = round(InFrameS_AA/TotalS_AA*100,1),
                     Strand = STRAND,
                     AA=z)
  
  INFO[is.na(INFO)] <- 0
  INFO2 <- INFO %>% filter(InFrameC>=10,InFrameC_perc>50,InFrameS_perc>30) #TotalC>10,InFrameC_perc>50,InFrameS_perc>50
  # INFO2 <- INFO %>% filter(TotalC>=10)
  print(nrow(INFO2))
  INFO2
}

AA2 <- Length$aa[-c(1)]

InFrameCountSiteLists <- list()
for(i in seq_along(AA2)){
  InFrameCountSiteLists[[i]]=InFrameCountSites(z=AA2[i])
}

sum(unlist(length(InFrameCountSiteLists))) #138
ALL <- sum(unlist(lapply(seq_len(length(InFrameCountSiteLists)),function(x) nrow(InFrameCountSiteLists[[x]]))))
ALL #7058

InFrameCountSiteDF <- do.call("rbind", InFrameCountSiteLists)
unique(InFrameCountSiteDF$AA) #Those peptide length identified with TuORFs

###Kallisto with 2 max isoforms output, provide the tpm values
Kallisto_2MaxIso_out <- read.delim(file="~/Desktop/CTRL_TPC/Kallisto01272023/2MaxIsoforms/abundance.tsv",header=T,stringsAsFactors=F,sep="\t")
head(Kallisto_2MaxIso_out,10)
Kallisto_2MaxIso_out$tx_id <- Kallisto_2MaxIso_out$target_id
Kallisto_2MaxIso_out$gene_id <- substr(Kallisto_2MaxIso_out$tx_id,1,9)
Kallisto_2MaxIso_out2 <- Kallisto_2MaxIso_out %>% dplyr::select(tx_id,tpm)

# Calculate the counts in different lengths of uORFs
TuORFs_df <- data.frame()
length(unique(InFrameCountSiteDF$AA))
for(i in unique(InFrameCountSiteDF$AA)){
  # print(i) #Last one is 335
  fiveUTR_ORFs_xaa <- fiveUTR_ORFs[uORF_length==(i+1)*3]
  TuORF <- InFrameCountSites(z=i)
  TuORF$tx_id <- sub("_.*","",TuORF$ORF_id)
  TuORF$chr <- substr(TuORF$ORF_id,3,3)
  
  # Define coordinates for minimum uORF start or stop on transcripts
  ORF_id <- TuORF$ORF_id
  Tx_id <- gsub("\\_.*","",ORF_id)
  
  uORF_Tx_START <- function(x){
    TxCord <- mapToTranscripts(fiveUTR_ORFs_xaa[[ORF_id[x]]], exonByTx[Tx_id[x]], ignore.strand = FALSE)
    start(TxCord[1])
  }
  
  uORF_Tx_END <- function(x){
    TxCord <- mapToTranscripts(fiveUTR_ORFs_xaa[[ORF_id[x]]], exonByTx[Tx_id[x]], ignore.strand = FALSE)
    end(TxCord[length(TxCord)])
  }
  
  #Find tx ranges for minimum uORFs
  uORF_Tx_START_data <- sapply(seq_len(length(ORF_id)),function(x) uORF_Tx_START(x))
  uORF_Tx_END_data <- sapply(seq_len(length(ORF_id)),function(x) uORF_Tx_END(x))
  TuORF$uORF_Tx_range <- paste0(uORF_Tx_START_data,"-",uORF_Tx_END_data)
  
  #Find genomic ranges for minimum uORFs
  uORF_Ge_START <- function(x){
    ifelse(as.character(strand(fiveUTR_ORFs[TuORF$ORF_id[x]])[[1]])[1]=="+",start(unlist(fiveUTR_ORFs[TuORF$ORF_id[x]]))[1],end(unlist(fiveUTR_ORFs[TuORF$ORF_id[x]]))[1])
  }
  
  uORF_Ge_END <- function(x){
    ifelse(as.character(strand(fiveUTR_ORFs[TuORF$ORF_id[x]])[[1]])[1]=="+",max(end(unlist(fiveUTR_ORFs[TuORF$ORF_id[x]]))),min(start(unlist(fiveUTR_ORFs[TuORF$ORF_id[x]]))))
  }
  
  uORF_Ge_START_data <- sapply(seq_len(length(ORF_id)),function(x) uORF_Ge_START(x))
  uORF_Ge_END_data <- sapply(seq_len(length(ORF_id)),function(x) uORF_Ge_END(x))
  
  TuORF$uORF_Ge_range <- paste0(uORF_Ge_START_data,"-",uORF_Ge_END_data)
  
  ###Remove redundant uORFs by genomic ranges
  TuORF2 <- left_join(TuORF,Kallisto_2MaxIso_out2,by="tx_id")
  TuORF2b <- TuORF2 %>% group_by(chr,uORF_Ge_range) %>% 
    mutate(tx_id_combine = toString(tx_id)) %>% 
    as.data.frame()
  
  ### Add peptide sequences
  fiveUTR_ORFs_seq_df <- fiveUTR_ORFs_seq %>% as.data.frame()
  fiveUTR_ORFs_peptide_df <- translate(fiveUTR_ORFs_seq) %>% as.data.frame()
  colnames(fiveUTR_ORFs_peptide_df) <- c("peptide_seq")
  fiveUTR_ORFs_peptide_df$ORF_id <- rownames(fiveUTR_ORFs_peptide_df)
  TuORF2c <- left_join(TuORF2b,fiveUTR_ORFs_peptide_df,by = "ORF_id")
  saveRDS(TuORF2c,paste0("~/Desktop/CTRL_TPC/Tiny_uORF_data/2Max_isoform_0.5TPM_Length_",i,"aa_ORFs.RData"))
  TuORFs_df<- rbind(TuORFs_df,TuORF2c)
}

TuORFs_df$gene_id <- substr(TuORFs_df$tx_id,1,9)
saveRDS(TuORFs_df,file="~/Desktop/CTRL_TPC/data/Rdata/uORF_gt2aa_Mar03-2023_10Ribo_50inFrame_30psite.RData")
TuORFs_df <- readRDS("~/Desktop/CTRL_TPC/data/Rdata/uORF_gt2aa_Mar03-2023_10Ribo_50inFrame_30psite.RData")

TuORFs_df2 <- TuORFs_df %>% 
  group_by(uORF_Ge_range,peptide_seq) %>% 
  dplyr::mutate(ORF_id_combine = toString(ORF_id)) %>% 
  slice_max(order_by = tpm, n = 1) %>% 
  as.data.frame()
nrow(TuORFs_df2) #6169

TuORFs_df2.5 <- TuORFs_df %>% 
  group_by(uORF_Ge_range) %>% 
  dplyr::mutate(ORF_id_combine = toString(ORF_id)) %>% 
  slice_max(order_by = tpm, n = 1) %>% 
  as.data.frame()
nrow(TuORFs_df2.5) #6163

#Find uORFs for same genome start and stop, but different AA sequences
setdiff(TuORFs_df2$ORF_id,TuORFs_df2.5$ORF_id)
# [1] "AT1G05160.2_4" "AT4G31350.2_3" "AT2G40830.3_1"
# [4] "AT5G53500.2_1" "AT1G07700.2_2" "AT5G08750.2_1"

TuORFs_df3 <- TuORFs_df %>% 
  group_by(uORF_Ge_range) %>% 
  mutate(ORF_id_combine = toString(ORF_id)) %>% 
  as.data.frame()

nrow(TuORFs_df3) #7058
write.xlsx(TuORFs_df2,file ="~/Desktop/CTRL_TPC/Tables/TuORFs_2max_isoforms_gl_10_inFrame_reads_50inFrame_30sites_6169_Mar3_2023.xlsx" )
write.xlsx(TuORFs_df3,file ="~/Desktop/CTRL_TPC/Tables/TuORFs_2max_isoforms_gl_10_inFrame_reads_50inFrame_30sites_7058_Mar3_2023_duplicated_uORF_not_removed.xlsx" )

```
