For the files required for this analysis see [here](https://data.mendeley.com/datasets/89j7snbm2r/draft?a=7b3a1001-2727-4a31-8ae5-784337dcc582)  
Use Figure 9A as an example.  

```
library(dplyr)
library(DataScienceR)
library(ggplot2)

ORF_max_filt <- read.delim(file="~/Desktop/CTRL_v1/ORFs_max_filt_expressed",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
ORF_max_filt2 <- ORF_max_filt %>% filter(!is.na(category)) %>% filter(category!="Overl_dORF")
pairwise_ks_test(value = ORF_max_filt2$ORF_Psit_pct_in_frame,group = ORF_max_filt2$category)
pIF <- ggplot(ORF_max_filt2, aes(x=ORF_Psit_pct_in_frame,color=category))+
  stat_ecdf(geom = "step")+
  xlim(0.5,0.999)+
  labs(tag = "Figure 9A") +
  theme_classic()
pIF
```
![image](https://github.com/hsinyenwu/ORFeome/assets/4383665/732812e9-65df-4ca7-a101-feac58a57b05)

Next calculate pairwise Kolmogorov–Smirnov test and obtain comtiple comparison letters.  
```
A=pairwise_ks_test(value = ORF_max_filt2$ORF_Psit_pct_in_frame,group = ORF_max_filt2$category)
A
library(multcompView)
multcompLetters(A*factorial(4-1),
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)
# ORFs_ccds    ncORFS      uORF      dORF 
#      "a"        "b"       "c"      "bc" 
```


