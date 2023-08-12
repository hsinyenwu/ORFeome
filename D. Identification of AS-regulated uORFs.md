Find isoform-specific uORFs.
```
#1: Find which gene has more than 1 isoform in the gtf
GTF_path <- "~/xxx/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf"
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF(file = GTF_path,format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByTx <- exonsBy(txdb, by="tx", use.names=T)
num_tx <- names(exonByTx)
num_tx_gene <- substr(num_tx,1,9)
df <- data.frame(num_tx,num_tx_gene)
df3 <- df %>% dplyr::count(num_tx_gene) 
table(df3$n)
# 1     2     3 
# 16513  2663     2 
df4 <- df3 %>% filter(n==2)
TuORFs_df2$nchar_tx_id_combine <- nchar(TuORFs_df2$tx_id_combine)
Iso_specific_TuORFs_df2 <- TuORFs_df2 %>% filter(nchar_tx_id_combine<17) #means from 1 isoform
nrow(Iso_specific_TuORFs_df2) #5281
# Select only the genes with more than one isoform
Iso_specific_TuORFs_df3 <- Iso_specific_TuORFs_df2 %>% filter(gene_id%in%df4$num_tx_gene)
nrow(Iso_specific_TuORFs_df3) #727
write.xlsx(Iso_specific_TuORFs_df3,file ="~/xxx/TuORFs_2max_isoforms_gl_10_inFrame_reads_50inFrame_30sites_727_isoform_specific.xlsx" )

#Load 5'UTR info
# GTF_path <- "~/xxx/Araport11_20181206_max2_expressed_isoform_TPM0.5_rm_diff_CDSstart.gtf"
# FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
# txdb <- makeTxDbFromGFF(file = GTF_path,format="gtf", dataSource="Araport11",organism="Arabidopsis")
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)

# Only retain both isoforms with  5'UTR ranges with they have different splice junctions 
# i.e., remove isoforms with identical 5'UTR splice junctions or no splice junctions 
Gene_Tx <- data.frame(gene_id=substr(names(fiveUTR),1,9),tx_id=names(fiveUTR))
Gene_Tx_list <- split(Gene_Tx,f=Gene_Tx$gene_id)
length(Gene_Tx_list) #18363
Gene_Tx_num <- sapply(1:length(Gene_Tx_list),function(x) nrow(Gene_Tx_list[[x]]))
Gene_Tx_list <- Gene_Tx_list[Gene_Tx_num==2] #both isoforms need to have 5'UTR to be considered further
length(Gene_Tx_list) #2571

# Determine if splice junction is the same
Splicing_same <- c()
for(i in 1:2571){
  if(i%%200==0) print(i)
  A <- gaps(unlist(fiveUTR[Gene_Tx_list[[i]][1,2]]),start=NA)
  B <- gaps(unlist(fiveUTR[Gene_Tx_list[[i]][2,2]]),start=NA)
  Splicing_same[i] <- identical(A,B)
}

# Remove those with same 
Gene_Tx_list2 <- Gene_Tx_list[Splicing_same==FALSE]
length(Gene_Tx_list2) #1140

Iso <- read.xlsx("~/xxx/TuORFs_2max_isoforms_gl_10_inFrame_reads_50inFrame_30sites_727_isoform_specific.xlsx")
Iso2 <- Iso %>% filter(gene_id %in% names(Gene_Tx_list2))
nrow(Iso2) #594
nrow(Iso) #727

library(stringr)
dfIso <- df %>% filter(num_tx_gene %in% Iso2$gene_id) %>% group_by(num_tx_gene) %>% dplyr::summarize(tx_isoforms = str_c(num_tx, collapse = ", ")) %>% as.data.frame()
head(dfIso)

Gene_description <- read.xlsx("~/xxx/TAIR_gene_description_w_gene_id_Dec2020.xlsx",sheet=1)
Iso2_ginfo <- dplyr::left_join(Iso2,Gene_description,by="gene_id")
head(Iso2_ginfo)
list_of_datasets <- list("A_2_isofroms" = dfIso, "B_Iso_TuORFs" = Iso2_ginfo)
write.xlsx(list_of_datasets,file="~/xxx/TableS7_594_isoform_specific_TuORFs.xlsx")

```
