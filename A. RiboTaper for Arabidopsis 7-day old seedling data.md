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

Step 5:   
Next, select the i, x, y, o, u, s class from the newly assembled gtf  
Here is just an example, you can select the classes according to analysis needs.
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
write.table(gtf_b,"/xxx/Araport11+CTRL_20181206.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
```
Step 6:
Map the RNA-seq data again as above. Run RSEM (quantify gene expression in gene and isoform levels) and only use isoforms with > 0 TPM.
```
# Here only shows the RSEM step
# With RSEM/1.3.1,STAR/2.6.0c
# Create RSEM index file

newINDEX=/xxx/RSEM_RNA_index
FASTA=/xxx/TAIR10_chr_all_2.fas #Watch out how the chromosome number is named! Here is 0, 1, 2, ...
GTF=/xxx/Araport11+CTRL_20181206.gtf

# Generate new RNA index with updated gtf for RSEM
echo "Generate new index with updated gtf for RSEM"

mkdir -p $newINDEX
cd $newINDEX

rsem-prepare-reference \
--gtf $GTF \
--star --star-path /opt/software/STAR/2.6.0c-foss-2018a/bin \
--star-sjdboverhang 99 \
-p 8 \
$FASTA \
RNA
```
Next, run RSEM
```
# Run RSEM
# DATA is the sample variable  
IPNUT_RNA=/xxx/STAR2
OUTNUT_RNA=/xxx/RSEM_Whole_Tx
INDEX_RNA=/xxx/RSEM_RNA_index

rsem-calculate-expression \
 --paired-end --bam --no-bam-output -p 10 --time \
 --strandedness reverse \
 --alignments $IPNUT_RNA/$DATA/"star_"$DATA"_"Aligned.toTranscriptome.out.bam $INDEX_RNA/RNA $DATA
```
Next use R to select expressed isoforms and make a new gtf file.
```
#Load RSEM output files
D1 <- read.delim("/xxx/D1.isoforms.results",header=T,sep="\t",stringsAsFactors = F, quote="")
#Also load D2, D3, D5,D6,D7 files 

TPM <- data.frame(gene_id=D1$gene_id, transcript_id=D1$transcript_id, D1_TPM=as.numeric(D1$TPM), D2_TPM=as.numeric(D2$TPM), D3_TPM=as.numeric(D3$TPM), D5_TPM=as.numeric(D5$TPM), D6_TPM=as.numeric(D6$TPM), D7_TPM=as.numeric(D7$TPM),stringsAsFactors=F)
TPM$TPM_mean = rowMeans(TPM[,3:8])

library(dplyr)
#Select genes with > 0 TPM
TPM2 = TPM %>% filter(TPM$TPM_mean>0)
nrow(TPM2) #46826 
length(unique(TPM2$gene_id)) #29439
#How many new transcripts are over TPM 0
sum(grepl("MSTRG",TPM2$gene_id)) #3314

gtf0 <-read.delim("/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/Araport11+CTRL_20181206.gtf",header=F,sep="\t",stringsAsFactors = F,quote = "")
gtf <- gtf0 %>% filter(V3 %in% c("mRNA","exon","CDS"))
gtf$tx_id <- unname(sapply(gtf$V9,function(x)  gsub("\"","",sub(pattern = "transcript_id ",replacement = "",grep("transcript_id",unlist(strsplit(x,";")),value=T)))))
gtf$tx_id <-gsub(" ", "", gtf$tx_id, fixed = TRUE)
gtf$V9[which(gtf$V3=="CDS")] <- paste0(gtf$V9[which(gtf$V3=="CDS")], " gene_biotype \"protein_coding\";")

gtf$biotype <- unname(sapply(gtf$V9,function(x)  gsub("\"","",sub(pattern = " gene_biotype ",replacement = "",grep(" gene_biotype ",unlist(strsplit(x,";")),value=T)))))
head(gtf)

nrow(gtf) #665559
gtf2 <- gtf %>% filter(tx_id %in% TPM2$transcript_id)
nrow(gtf2) #572270
gtf2$gene_id <- sub("^(.*)[.].*", "\\1", gtf2$tx_id)
head(gtf2,60)
gtf2_s <- split(gtf2,gtf2$gene_id)

gene_V9 <- function(x) {
        gene_row=x[which(x$V3=="mRNA")[1],]
        gene_row$V3="gene"
        gene_row$V4=min(x[which(x$V3=="mRNA"),4])
        gene_row$V5=max(x[which(x$V3=="mRNA"),5])
        gene_row$V9=paste0("gene_id \"", x$gene_id[1],"\"; ","gene_biotype \"",x$biotype[1],"\";")
        rbind(gene_row,x)}


gtf2_s2 <- lapply(gtf2_s, function(x) gene_V9(x))
head(gtf2_s2,2)
tail(gtf2_s2,2)
gtf2_s3 <- bind_rows(gtf2_s2)
nrow(gtf2_s3) #601709
#There was no space before "transcript_biotype",fix it
gtf2_s3$V9 <- gsub("transcript_biotype"," transcript_biotype",gtf2_s3$V9)
tail(gtf2_s3,60)
#Araport11+CTRL_20181206.gtf is for the isoform TPM > 0.5
#write.table(gtf2_s3[,1:9],"/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/Araport11+CTRL_20181206_TPM0.5.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
#Araport11+CTRL_20181206_TPM1.gtf is for TPM > 1
write.table(gtf2_s3[,1:9],"/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/Araport11+CTRL_20181206_expressed.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
```
Step 7.
Run STAR again with the new annotation for RNA-seq.  
Run STAR (1st time) with the new annotation for Ribo-seq.  

```
# With STAR/2.6.0c
# STAR index for Ribo-seq
# Path to folder for the output: newINDEXRibo
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $newINDEXRibo \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 34 \

# STAR mapping for Ribo-seq
STAR --runThreadN 10 \
     --genomeDir $newINDEXRibo \
     --alignEndsType EndToEnd \
     --readFilesCommand zcat \
     --readFilesIn $INPUT/$DATA.noContam4.fastq.gz \
     --alignIntronMax 5000 \
     --alignIntronMin 15 \
     --outFilterMismatchNmax 2 \
     --outFilterMultimapNmax 20 \
     --outFilterType BySJout \
     --alignSJoverhangMin 4 \
     --alignSJDBoverhangMin 1 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outSAMmultNmax 1 \
     --outMultimapperOrder Random \
     --outFileNamePrefix "star_"$DATA"_" \
```
```
# STAR mapping for RNA-seq
# RNA-seq index (newINDEXRNA) can be created by above code and the new gtf
STAR --runThreadN 10 \
     --genomeDir $newINDEXRNA \
     --readFilesCommand zcat \
     --readFilesIn $INPUT/$DATA.r1.fastq.gz $INPUT/$DATA.r2.fastq.gz \
     --alignIntronMax 5000 \
     --alignIntronMin 15 \
     --outFilterMismatchNmax 2 \
     --outFilterMultimapNmax 20 \
     --outFilterType BySJout \
     --alignSJoverhangMin 2 \
     --alignSJDBoverhangMin 1 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outSAMmultNmax 1 \
     --outMultimapperOrder Random \
     --outFileNamePrefix "star_"$DATA"_" \

```
Step 8.  
Make RiboTaper annotation  
```
#With R/3.4.3-X11-20171023, RiboTaper_v1.3, bedtools code from RiboTaper website
BEDTOOL=/xxx/bedtools_dir

# The code should look like:
# create_annotation_files.bash <gencode_gtf_file> <genome_fasta_file(indexed)> <use_ccdsid?> <use_appris?> <dest_folder> <bedtools_path> <scripts_dir>

$xxx/RiboTaper_v1.3/scripts/create_annotations_files.bash $GTF $FASTA false false $OUTPUT $BEDTOOL $xxx/RiboTaper_v1.3/scripts/
```
Step 9.  
Merge output and run RiboTaper  
```
#SAMtools/1.8
# next merge (with samtools merge) the RNA-seq bam files to one bam file and merge Ribo-seq bam files to another bam file. Sort those files to genome coordinates (with samtools sort) and create the index (with samtools index)
simplified example code:
#Merge D1, D2, D3...files
samtools merge -f $INPUT_ribo/ribo_CTRL_merged.bam $INPUT_ribo/D1/star*sortedByCoord.out.bam $INPUT_ribo/D2/star*sortedByCoord.out.bam ...
#Sort by genome coordinates
samtools sort -o $INPUT_ribo/ribo_CTRL_merged2.bam $INPUT_ribo/ribo_CTRL_merged.bam
#Change name
mv $INPUT_ribo/ribo_CTRL_merged2.bam $INPUT_ribo/ribo_CTRL_merged.bam
#Create index
samtools index $INPUT_ribo/ribo_CTRL_merged.bam
```
Run RiboTaper
```
TAPER: Path to RiboTaper code e.g., /xxx/Software/RiboTaper_v1.3/scripts
RNA, RIBO: Path to the merged RNA and Ribo-seq folder
ANNO: Path to the RiboTaper annotation folder (created above)
BED: Path to the bedtools

#RiboTaper Usage: ./Ribotaper.sh <Ribo_bamfile> <RNA_bamfile> <annotation_dir> <comma-sep_read_lenghts_ribo> <comma-sep_cutoffs> <scripts_dir> <bedtools_dir>  <n_cores> 

$TAPER/Ribotaper.sh $RIBO/ribo_CTRL_merged.bam $RNA/RNA_CTRL_merged.bam $ANNO 24,25,26,27,28 8,9,10,11,12 $TAPER $BED 8
```
Step 10. Make gtf with only CDS ranges and quantify with RSEM. 
Use code similar to above R code and RSEM index/quant code. 


