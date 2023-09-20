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
# Run stringtie to assemble transcripts
echo "### Run stringtie to assemble transcripts for each sample"
# With StringTie/1.3.5

cd $assembledGTF
stringtie --rf -p 16 -G $GTF -o $DATA.gtf -l $DATA $OUTPUT2/$DATA/star*sortedByCoord.out.bam
```



