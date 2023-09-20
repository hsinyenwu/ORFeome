Step 1: Create Bowtie2 index for removing contamination sequences (rRNA, tRNA, snRNA and snoRNA  
```
#With Bowtie2/2.3.4
#CF is the path to the contamination sequence folder
#Output is named Contam4
bowtie2-build $CF/Araport11_201606_rRNA.fasta,$CF/Araport11_201606_tRNA.fasta,$CF/Araport11_201606_snRNA.fasta,$CF/Araport11_201606_snoRNA.fasta Contam4
```
Step 2: preprocessing Ribo-seq reads  
```
#With Bowtie2/2.3.4, SAMtools/1.8, BEDTools/2.26.0
#For each Ribo-seq fasta file, run the following. You can also create a loop to run these, or submit multiple jobs with a cluster.
#Folder path variables: INPUT, OUTPUT, ContamIndex (from the Contam4)
#Sample variable (the name of each sample, : DATA

#*********************#
#direct to output directory
cd $OUTPUT
#remove rRNA, tRNA, snRNA/snoRNA & repeat sequences
echo "### rRNA, tRNA, snRNA, snoRNA & repeat sequences ###"
echo "  $DATA"
bowtie2 -L 20 -p 8 -x $ContamIndex -1 $INPUT/$DATA"_r1.fastq.gz" -2 $INPUT/$DATA"_r2.fastq.gz" -S $DATA.mapped_and_unmapped.sam
#-L 20 means sequence seed lindex length is 20

echo "sam to bam"
samtools view -bS -o $DATA.mapped_and_unmapped.bam $DATA.mapped_and_unmapped.sam

#Extract unmapped reads (This is what we want)
# -F 256 is discarding "secondary alignments"
# -f 12 extract unmapped reads whose mates are also unmapped
echo "### Extract unmapped reads ###"
echo "  $DATA"
samtools view -b -f 12 -F 256 -o $DATA.bothEndsUnmapped.bam $DATA.mapped_and_unmapped.bam

#Split paired-end reads into separated fastq files .._r1 .._r2
echo "### Split paired-end reads into separated fastq files ###"
samtools sort -n -o $DATA.bothEndsUnmapped_sorted.bam $DATA.bothEndsUnmapped.bam

bamToFastq -i $DATA.bothEndsUnmapped_sorted.bam -fq $DATA.r1.fastq -fq2 $DATA.r2.fastq

gzip -f $DATA.r1.fastq
gzip -f $DATA.r2.fastq
```
Next do the same for the RNA-seq files. 

Step 3. First STAR mapping for RNA-seq data (this will be used for reference-guided transcriptome assembly)  
Create index:  
```
#With STAR/2.6.0c
#Define path to FASTA and GTF file: FASTA, GTF
#Path to the index: newINDEX
mkdir -p $newINDEX
cd $newINDEX

echo "Generate index with  for STAR"
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $newINDEX \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 99 \
```
Map RNA-seq reads:  
```
# With STAR/2.6.0c, SAMtools/1.8
# Path to folders: INPUT, OUTPUT, OUTPUT2
# From previous step: starIndex1
# Path for splice junction file: SJ

echo "### Run STAR to generate genome-mapped bam files for stringtie ###"
echo "see $assembledGTF for output gtf"
cd $OUTPUT/$DATA

#1st pass STAR alignment        
STAR --runThreadN 10 \
--genomeDir $starIndex1 \
--readFilesCommand zcat \
--readFilesIn $INPUT/$DATA.r1.fastq.gz $INPUT/$DATA.r2.fastq.gz \
--alignIntronMax 5000 \
--alignIntronMin 15 \
--outFilterMismatchNmax 2 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 2 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFileNamePrefix "star_"$DATA"_" \

rm "star_"$DATA"_"Aligned.toTranscriptome.out.bam

cd $OUTPUT2/$DATA

#2nd pass STAR alignment
STAR --runThreadN 10 \
--genomeDir $starIndex1 \
--sjdbFileChrStartEnd $SJ \
--readFilesCommand zcat \
--readFilesIn $INPUT/$DATA.r1.fastq.gz $INPUT/$DATA.r2.fastq.gz \
--alignIntronMax 5000 \
--alignIntronMin 15 \
--outFilterMismatchNmax 2 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 2 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFileNamePrefix "star_"$DATA"_" \

samtools index "star_"$DATA"_"Aligned.sortedByCoord.out.bam
rm "star_"$DATA"_"Aligned.toTranscriptome.out.bam
```

Step 4:
StringTie to assemble transcripts for each file.  
```
# With StringTie/1.3.5
# Run stringtie to assemble transcripts for each sample
echo "### Run stringtie to assemble transcripts for each sample"

cd $assembledGTF
stringtie --rf -p 16 -G $GTF -o $DATA.gtf -l $DATA $OUTPUT2/$DATA/star*sortedByCoord.out.bam
```
Merge transcripts from all samples:
```
echo "### Merge transcripts from all samples ###"
stringtie --merge -p 12 -T 0.05 -G $GTF -o CTRL_Araport11_merged.gtf $assembledGTF/mergeList.txt

echo "### Compare transcripts with the annoation ###"
# With gffcompare
# xxx path to the gffcompare code
/xxx/gffcompare-0.11.2/gffcompare -V -r $GTF -o CTRL_Araport11 CTRL_Araport11_merged.gtf #this creates CTRL_Araport11.annotated.gtf
```
Nest, select the i, x, y, o, u, s class from the newly assembled gtf
```
Last login: Wed Sep 20 09:19:52 on ttys000

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
(base) Wus-MacBook-Pro-3:Proteomics wu$ ssh -XY larrywu@hpcc.msu.edu
Last login: Fri Aug 25 10:36:39 2023 from 149.102.252.46
S                 __ _____ _____ ____      __  __ ____  _____ _____              
                /_// ___// ___// _  ;    / /_/ // _  ;/ ___// ___/              
               / // /__ / _/_ / _,.'    / __  // _,.'/ /__ / /__                
              /_//____//____//_//_!    /_/ /_//_/   /____//____/                
________________________________________________________________________________
                                                                                
          Welcome to Michigan State's High Performance Computing Center         
                    ** Unauthorized access is prohibited **                     
________________________________________________________________________________
                                                                                
    For GPU development please use green nodes.                                 
                                                                                
          Development Nodes (usage)                                             
      --------------------------------                                          
S    dev-amd20-v100 (low) dev-amd20 (low)                             
    dev-intel14-k20 (low) dev-intel14 (low)                                  
    dev-intel16-k80 (low) dev-intel16 (low)                                                               
    dev-intel18 (low)                                                                                                  
                                                                                
      |--------------------------------------------------------------------------------------| 
      |                                                                                      | 
      |   New blog and documenation sites:                                                   | 
      |     https://docs.icer.msu.edu                                                        | 
      |     https://blog.icer.msu.edu                                                        | 
      |                                                                                      | 
      |--------------------------------------------------------------------------------------| 
                                                                               
ssh dev-amd20________________________________________________________________________________
[larrywu@gateway-03 ~]$ ssh dev-amd20
Last login: Wed Sep 20 09:20:06 2023 from gateway-01.dmz

===
Please note that processes on development nodes are limited to two hours of 
CPU time; for longer-running jobs, please submit to the queue.

Development nodes are a shared system; for information about performance 
considerations please see: https://docs.icer.msu.edu/development_nodes/
===

[larrywu@dev-amd20 ~]$ cd ABA
[larrywu@dev-amd20 ABA]$ cd code/
[larrywu@dev-amd20 code]$ ll
total 1005
-rw-r----- 1 larrywu biochem     740 Dec  6  2018 1.0_bowtie2_index.sb
-rwxr-x--- 1 larrywu biochem     361 Dec  6  2018 1.1a_RNA_preprocessing_submission.sh
-rw-r----- 1 larrywu biochem    1886 Dec  6  2018 1.1_RNA_preprocessing.sb
-rwxr-x--- 1 larrywu biochem     355 Dec  6  2018 1.2a_ribo_preprocessing_submission.sh
-rw-r----- 1 larrywu biochem    1841 Dec  6  2018 1.2_ribo_preprocessing.sb
-rw-r----- 1 larrywu biochem    1101 Dec 11  2018 1.5_STAR_index_Araport.sb
-rw-r----- 1 larrywu biochem    2167 Dec 11  2018 2.0_RNA_mapping.sb
-rwxr-x--- 1 larrywu biochem     645 Dec  7  2018 2.0_RNA_mapping_submission.sh
-rw-r----- 1 larrywu biochem    1100 Dec 11  2018 2.1_gtf_merge_compare.sb
-rw-r----- 1 larrywu biochem     565 Dec 11  2018 3.0a_Rscript_gtf_processing.sb
-rw-r----- 1 larrywu biochem    3724 Dec 11  2018 3.0_stringtie_merged_gtf_analysis.R
-rw-r----- 1 larrywu biochem    1313 Dec 11  2018 3.1a_Index_Kallisto_RNAseq.sb
-rw-r----- 1 larrywu biochem    1013 Dec 11  2018 3.1b_Quant_Kallisto_RNAseq.sb
-rwxr-x--- 1 larrywu biochem     292 Dec 11  2018 3.1b_Quant_Kallisto_submission.sh
-rw-r----- 1 larrywu biochem    1149 Dec 11  2018 3.1k_RSEM_RNA_index.sb
-rw-r----- 1 larrywu biochem    1033 Dec 11  2018 3.1m_RSEM_RNA_quant.sb
-rwxr-x--- 1 larrywu biochem     285 Dec 11  2018 3.1m_RSEM_RNA_quant_submission.sh
-rw-r----- 1 larrywu biochem    4415 Dec 12  2018 3.1n_get_highly_expressed_isoforms_and_make_gtf.R
-rw-r----- 1 larrywu biochem    1326 Dec 12  2018 3.2a_ribotaper_annotation.sb
-rw-r----- 1 larrywu biochem    1562 Dec 12  2018 3.3a_STAR_index_Araport.sb
-rw-r----- 1 larrywu biochem    1544 Dec 12  2018 3.3b_ribo_mapping.sb
-rw-r----- 1 larrywu biochem    1527 Dec 12  2018 3.3b_RNA_mapping.sb
-rwxr-x--- 1 larrywu biochem     865 Dec 12  2018 3.3b_RNA_ribo_mapping_submission.sh
-rwxr-x--- 1 larrywu biochem     265 Dec 12  2018 3.3c_merge_submission.sh
-rw-r----- 1 larrywu biochem    1164 Dec 12  2018 3.3c_ribo_ABA60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1329 Dec 12  2018 3.3c_ribo_ABA60+DMSO60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1171 Dec 12  2018 3.3c_ribo_DMSO60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1138 Dec 12  2018 3.3c_RNA_ABA60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1308 Dec 12  2018 3.3c_RNA_ABA60+DMSO60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1215 Dec 12  2018 3.3c_RNA_DMSO60_bam_merge.sb
-rw-r----- 1 larrywu biochem    1420 Dec 13  2018 3.4_ABA_RiboTaper_template_20181206.sb
-rw-r----- 1 larrywu biochem    1328 Dec 12  2018 3.4b_ABA60+DMSO60_RiboTaper_20181206.sb
-rw-r----- 1 larrywu biochem    1528 Dec 14  2018 3.4c_ABA60+DMSO60_RiboTaper_old_node_20181206.sb
-rwxr-x--- 1 larrywu biochem     351 Dec 12  2018 3.4_RiboTaper_submission_20181206.sh
-rw-r----- 1 larrywu biochem     116 Dec  7  2018 3.6_command_line_merge_P_sites_all.txt
-rw-r----- 1 larrywu biochem    2477 Dec  7  2018 4.0b_Index+Quantfy_Kallisto_tx.sb
-rw-r----- 1 larrywu biochem    1141 Nov 26  2018 5.1c_Tomato_ribo_CDS_quant.sb
-rw-r----- 1 larrywu biochem    2462 Dec 26  2018 6.0_Index+Quant_Ctrl+ABA_CDS_Kallisto.sb
-rw-r----- 1 larrywu biochem    1548 Dec 11  2018 old_3.1_Index+Quantfy_Kallisto_tx.sb
drwxr-x--- 2 larrywu biochem   16384 Dec 14  2020 out
-rw-r----- 1 larrywu biochem     267 Dec 12  2018 slurm-3713058.out
-rw-r----- 1 larrywu biochem     267 Dec 12  2018 slurm-3713089.out
-rw-r----- 1 larrywu biochem   28365 Dec 13  2018 slurm-3725903.out
-rw-r----- 1 larrywu biochem   28365 Dec 13  2018 slurm-4173534.out
-rw-r----- 1 larrywu biochem   28371 Dec 14  2018 slurm-4280714.out
-rw-r----- 1 larrywu biochem    5220 Dec 14  2018 slurm-4329883.out
-rw-r----- 1 larrywu biochem   45044 Dec 26  2018 slurm-5607405.out
-rw-r----- 1 larrywu biochem   32100 Dec 26  2018 slurm-5608694.out
-rw-r----- 1 larrywu biochem    1903 Dec  7  2018 Test_1.3_Index+Quantfy_Kallisto_tx.sb
-rw-r----- 1 larrywu biochem 2757472 Jan 31  2019 test.txt
drwxr-x--- 2 larrywu biochem    8192 Dec 12  2018 unused_code
[larrywu@dev-amd20 code]$ vi 1.0_bowtie2_index.sb
[larrywu@dev-amd20 code]$ vi 1.0_bowtie2_index.sb
[larrywu@dev-amd20 code]$ cd /mnt/home/larrywu/CTRL_arabidopsis/code_TPM0.25/3.0_stringtie_merged_gtf_analysis.R
-bash: cd: /mnt/home/larrywu/CTRL_arabidopsis/code_TPM0.25/3.0_stringtie_merged_gtf_analysis.R: No such file or directory
[larrywu@dev-amd20 code]$ vi 3.0_stringtie_merged_gtf_analysis.R

  gene_V4 <- min(mRNA$V4)
  gene_V5 <- max(mRNA$V5)
  gene <- mRNA[1,]
  gene$V3 <-  "gene"
  gene$V4 <- gene_V4
  gene$V5 <- gene_V5
  gene$V9 <- paste0("gene_id \"",x$gene_id[1],"\"; gene_biotype \"ncRNA\"; class_code ",x$class_code[1],";")
  rbind(gene,x)[,1:9]
}

# load dataframe from known genes
Araport11_gtf <- read.delim("/xxx/Araport11_20181206.gtf",header=F,sep="\t",stringsAsFactors = F,quote="",skip=0)

# data.frame for novel genes
MS_gtf2 <- lapply(MS_gtf_list,function(x) gene_rows(x))

# Dataframe from MS_gtf
MS_new_gtf <- do.call("rbind", MS_gtf2)
nrow(MS_new_gtf) #17596
MS_new_gtf[which(MS_new_gtf$V3=="transcript"),3] <- "mRNA"

# Combine new and original gtfs
gtf_b <- rbind(Araport11_gtf,MS_new_gtf)

# Save the gtf file
write.table(gtf_b,"/xxx/Araport11+CTRL_20181206.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
```
Step 5: 
```
library(dplyr)
rm(list=ls())
#load the gtf from stringtie output
#ABA_Araport11_20181206.annotated.gtf
gtf <- read.delim("/xxx/Araport11+CTRL_20181206.gtf",header=F,sep="\t",stringsAsFactors = F,quote="",skip=0)

#split by transcript
i1 <- gtf$V3 == "transcript"
gtfL <- split(gtf, cumsum(i1))
#extract class_code
gtfL2 <- sapply(gtfL, function(x) gsub(pattern = " class_code ", replacement = "", grep("class_code",unlist(strsplit(x[1,9],";")),value = T)))
#select noval class_code i, x, y, o, u, s
gtfL3 <- c(gtfL[unname(which(gtfL2=="\"i\"" ))],gtfL[unname(which(gtfL2=="\"o\"" ))],gtfL[unname(which(gtfL2=="\"x\"" ))],gtfL[unname(which(gtfL2=="\"y\"" ))],gtfL[unname(which(gtfL2=="\"u\"" ))],gtfL[unname(which(gtfL2=="\"s\"" ))])
head(gtfL3)
#order the list by names
gtfL4 <- gtfL3[order(as.numeric(names(gtfL3)))]
#make the list to data.frame
gtf_c <- bind_rows(gtfL4)

##############
# Prepare GTF for RiboTaper, STAR and RSEM
# Use stringtie assemble and merged gtf
# RiboTaper needs gene_biotype, which is lacking in the stringtie ouput gtf
# Need to add gene_biotype back for coding genes, for new genes, I will add ; gene_biotype "ncRNA";
# I did this on my mac. but could be done in HPCC
# Using quote="" in read.table is the key since gene_id and transcript_id need to be surrounded by double quotes

#extract transcript_id
gtf_c$transcript_id <- sapply(1:nrow(gtf_c), function(x) gsub("\"","", gsub("transcript_id ", replacement = "",unlist(strsplit(gtf_c$V9[x],"; "))[1])))
head(gtf_c$transcript_id)
#extract gene_id
gtf_c$gene_id <- sapply(1:nrow(gtf_c), function(x) gsub("\"","", gsub("gene_id ", replacement = "",unlist(strsplit(gtf_c$V9[x],"; "))[2])))
head(gtf_c$gene_id,20)
#extract class_code
gtf_c$class_code <- sapply(gtf_c[,9], function(x) gsub("\"","",gsub(pattern = " class_code ", replacement = "", grep("class_code",unlist(strsplit(x,";")),value = T))))
gtf_c$class_code <- unname(unlist(lapply(gtf_c$class_code, function(x) if(identical(x, character(0))) NA_character_ else x)))
table(gtf_c$class_code)

#   i    o    s    u    x    y 
# 193  188  121 1835 3300    4

for(i in 1:length(gtf_c$class_code)) {
  if(is.na(gtf_c$class_code[i])) {
    gtf_c$class_code[i] <- gtf_c$class_code[i-1]
  }
}

table(gtf_c$class_code)
#  i    o    s    u    x    y 
# 389  594  420 3891 6766   1

# Make each new gene to a list 
MS_gtf_list <- split(gtf_c,gtf_c$gene_id)

# A function to make the "gene" rows for new genes
gene_rows <- function(x) {
  x$V9 <- paste0(x$V9,"transcript_biotype \"ncRNA\"","; gene_biotype \"ncRNA\";")
  mRNA <- x[x$V3=="transcript",]
  gene_V4 <- min(mRNA$V4)
  gene_V5 <- max(mRNA$V5)
  gene <- mRNA[1,]
  gene$V3 <-  "gene"
  gene$V4 <- gene_V4
  gene$V5 <- gene_V5
  gene$V9 <- paste0("gene_id \"",x$gene_id[1],"\"; gene_biotype \"ncRNA\"; class_code ",x$class_code[1],";")
  rbind(gene,x)[,1:9]
}

# load dataframe from known genes
Araport11_gtf <- read.delim("/xxx/Araport11_20181206.gtf",header=F,sep="\t",stringsAsFactors = F,quote="",skip=0)

# data.frame for novel genes
MS_gtf2 <- lapply(MS_gtf_list,function(x) gene_rows(x))

# Dataframe from MS_gtf
MS_new_gtf <- do.call("rbind", MS_gtf2)
nrow(MS_new_gtf) #17596
MS_new_gtf[which(MS_new_gtf$V3=="transcript"),3] <- "mRNA"

# Combine new and original gtfs
gtf_b <- rbind(Araport11_gtf,MS_new_gtf)

# Save the gtf file
write.table(gtf_b,"/xxx/Araport11+CTRL_ixyous_20181206.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
```
  

