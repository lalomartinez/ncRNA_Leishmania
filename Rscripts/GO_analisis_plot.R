library(ggplot2)
library(ggpubr)




data1<-read.delim("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LdonBPK282_tritryp/DE_comparison_AMA_PRO_LdonBPK282_Tritryp_table.csv", 
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