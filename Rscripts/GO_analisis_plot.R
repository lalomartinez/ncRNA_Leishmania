#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

# to run this script you need as a input file: output file of GO terms enrichment

data1<-read.delim(args[1], 
                 sep = '\t', header = T)

data1$log.pvalue= -log(data1$p.value)
df1<-data.frame(data1$Description, data1$log.pvalue,data1$x)
colnames(df1)<-c("Description","Log.pvalue","X")

df1<-within(df1, 
            Description<- factor(Description, 
                                 levels=names(sort(table(Description), 
                                                   decreasing=TRUE))))

theme_set(
  theme_pubr() +
    theme(legend.position = "right")
)



plot1<- ggplot(df1, aes(x= reorder(Description, -Log.pvalue),Log.pvalue))+
  geom_bar(stat="identity", fill= '#4F779F', color='black',  position=position_dodge())+
  coord_flip()+
  labs(y="Log(p-value)", x="")+
  theme(axis.text.y = element_text(size=8))
plot1