library(goseq)
library(GO.db)
library(ggplot2)
library(ggpubr)
library(dplyr)
#read GO annotation file
gaf_ann= read.delim("~/Desktop/Leishmania_info_06-01-20/seq_Leishmania_Companion_GeneID/AA_Tritryp/gaf_files_Release35/TriTrypDB-35_LmajorFriedlin_GO.gaf", 
                    header = F, sep="\t")

final_ann<-data.frame(gaf_ann$V2, gaf_ann$V5)
final_ann<-final_ann[-1,]

# read Gene Length file and convert to numeric vector
gene_length= read.delim("~/Desktop/Leishmania_info_06-01-20/seq_Leishmania_Companion_GeneID/CDS_Tritryp/GC_content/TriTrypDB-35_LmajorFriedlin_AnnotatedCDSs.fasta.txt",
                        header = T, sep = "\t")
final_length<-data.frame(gene_length$ID,gene_length$Total.Count)
#final_length$gene_length.ID<-gsub("\\.1","", final_length$gene_length.ID)
fl<-as.numeric(as.vector(final_length$gene_length.Total.Count))
names(fl)<-final_length$gene_length.ID


# make vector of all genes
assayed_genes<- data.frame(final_length$gene_length.ID)
assayed_vector<-c(t(assayed_genes))


# make vector of DE genes or subset of interest genes
DE_genes<- read.delim("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/modules_comparison/LmajFriedlin_turquoise_GO.list", header = F)
DE_genes$V1<-gsub("\\.1","", DE_genes$V1)
DEG_vector<-c(t(DE_genes))
#Construct a new vector that adds a 0 next to every gene that is not in our DEG list and a 1 next to every gene that is in our DEG list.
gene.vector=as.integer(assayed_vector%in%DEG_vector)
names(gene.vector)=assayed_vector 

pwf=nullp(gene.vector, bias.data = fl)
head(pwf)

go_pwf<- goseq(pwf, gene2cat = final_ann, use_genes_without_cat=F)
head(go_pwf)


enriched_subset<-subset(go_pwf, go_pwf$over_represented_pvalue<.01)

# obtain a subset of enriched terms
enriched_GO=go_pwf$category[go_pwf$over_represented_pvalue<.01]
length(enriched_GO)



enriched_subset$log_pvalue= -log(enriched_subset$over_represented_pvalue)

df1<-data.frame(enriched_subset$term,enriched_subset$log_pvalue, enriched_subset$numDEInCat)
colnames(df1)<-c("Description","Log.pvalue","X")
df1$Description<-as.character(df1$Description)
df1$Description[is.na(df1$Description)]="Undefined"



df1<-within(df1, 
            Description<- factor(Description, 
                                 levels=names(sort(table(Description), 
                                                   decreasing=TRUE))))

theme_set(
  theme_pubr() +
    theme(legend.position = "right")
)


plot1<- ggplot(df1, aes(x= reorder(Description, -Log.pvalue),Log.pvalue))+
  geom_bar(stat="identity", fill= "#63BDD3", color='black',  position=position_dodge())+
  coord_flip()+
  labs(y="Log(p-value)", x="")+
  theme(axis.text.y = element_text(size=8))
plot1


go_pwf %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")



#obtain a file with all information about your enriched terms

capture.output(for(go in enriched_GO) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/modules_comparison/-vs-Ldon_turquoise-vs-LmajFriedlin_turquoise_sharedGO2Def.list")







#colores
#azul
#4F779F
# amarillo
#FFF759
#turquesa
#63BDD3
