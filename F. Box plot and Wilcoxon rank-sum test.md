

```
#######################
# DLA boxplot
#######################
rm(list=ls())
library(dplyr)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(grid)
library(EnvStats)
library(ggeasy)
#####################################
#Sample 1
# AT1G56300 
DLA <- read.xlsx("~/xxx/DLA_selected_for_plotting.xlsx",sheet=1)
colnames(DLA) <- c("Lines","DLA")
#test equal variance
# default p-value adjustement method: BY
stat.test <- DLA %>%
  wilcox_test(DLA ~ Lines) %>% 
  add_significance("p")
#  adjust_pvalue(method = "BY") %>%
stat.test

# Create a box plot
bxp <- ggboxplot(DLA, 
                 x = "Lines", y = "DLA", 
                 color = "Lines", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(0,3.5),
                 legend = "right") + 
  rotate_x_text(angle = 45) + 
  labs(tag = "") + 
  ylab("DL ratio")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  ggtitle("AT1G56300 \n DnaJ")+ geom_dotplot(binaxis='y', stackdir='centerwhole', dotsize=0.5) +
  ggeasy::easy_center_title() + theme(plot.title = element_text(face = "italic"))

yvalue <- function(ystart,stepsize,num) {
  Ev <- c()
  for(i in 1:num){
    Ev[i] <- ystart + stepsize*(i-1)
  }
  return(Ev)
}

vsteps <- yvalue(ystart=3,stepsize=0.25,num=1)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Lines", dodge = 0.8)

pDLA_AT1G56300 <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.015, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=0.1,size = 3.5)
pDLA_AT1G56300
```
