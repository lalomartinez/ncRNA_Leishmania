library(WGCNA)
options(stringsAsFactors = FALSE)

cts1 <- read.csv("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/normalizado/LbraM2903_TritrypVSncRNA_normLOG2.htseq",
                header = T,
                sep = "\t")

datExpr1 = as.data.frame(t(cts1[, -1]))
names(datExpr1) = cts1$gene_id
rownames(datExpr1) = names(cts1)[-1]

#Checking data for excessive missing values and identification of outlier samples

gsg= goodSamplesGenes(datExpr1,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)printFlush(paste("Removing genes:", paste(names(datExpr1)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]}

sampleTree1 = hclust(dist(datExpr1), method = "average")
sizeGrWindow(12,9)
plot(sampleTree1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 160, col = "red")

#charge aditional data

traitData1 = read.csv("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/experiments/metadata_LbraM2903_2.csv", sep = ',');
dim(traitData1)
names(traitData1)


# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr1);
traitRows1 = match(Samples, traitData1$SampleName);
datTraits1 = traitData1[traitRows1, -1];
rownames(datTraits1) = traitData1[traitRows1, 1];

collectGarbage();



#convert ordinal data to numeric data
#trait_factor= factor(traitData1$Class)
#unclass(trait_factor)
#traitData1$DevStage=as.numeric(trait_factor)

sampleTree3 = hclust(dist(datExpr1), method = "average")


traitColors1 = numbers2colors(datTraits1, signed = FALSE, colors = c("#FFFFFF", "#AE2B47"));
plotDendroAndColors(sampleTree3, traitColors1,groupLabels = names(datTraits1),main = "Sample dendrogram and development stage heatmap")

#save data for next steps 
save(datExpr1, datTraits1, file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/normalizado/LmajFriedlin_TritrypVSncRNA_normLOG2.htseq")

lnames1<- load(file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/LbraWGCNA_dataInput.RData")

powers1 = c(c(1:10), seq(from = 12, to=50, by=2))

sft1 = pickSoftThreshold(datExpr1, powerVector = powers1, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));text(sft1$fitIndices[,1], -sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],
                                              labels=powers1,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")


# Mean connectivity as a function of the soft-thresholding power
plot(sft1$fitIndices[,1], sft1$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft1$fitIndices[,1], sft1$fitIndices[,5], labels=powers1, cex=cex1,col="red")

#One-step network construction and module detection

net = blockwiseModules(datExpr1, power =12,TOMType = "unsigned", 
                       maxBlockSize = 150000,
                       minModuleSize = 30,reassignThreshold = 0, 
                       mergeCutHeight = 0.25,numericLabels = TRUE,
                       pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "LbraM2903_TOM",verbose = 3)
# contains the module assignment
net$colors

#contains the module eigengenes of the modules
net$MEs

#how many modules were identified and what the module sizes
#labeled 1 through N in order of descending size, with sizes ranging from X to Z genes. The label 0 is reserved for genes outside of all modules.
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#save module data for next steps
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

tail(names(datExpr1)[moduleColors=="turquoise"])




save(MEs1, moduleLabels1, moduleColors1, geneTree1,file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/LbraWGCNA_2modules_dataInput.RData")


# Load the expression and trait data saved in the first part
lnames3 = load(file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/LbraWGCNA_dataInput.RData");
#The variable lnames2 contains the names of loaded variables.
lnames3
# Load network data saved in the second part.
lnames3 = load(file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/LbraWGCNA_2modules_dataInput.RData");
lnames3

# Define numbers of genes andsamples
nGenes1 = ncol(datExpr1)
nSamples1 = nrow(datExpr1)

# Recalculate MEs with color labels
MEs01 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs1 = orderMEs(MEs01)
moduleTraitCor1 = cor(MEs1, datTraits1, use = "p")
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples1)

sizeGrWindow(10,9)
# Will display correlations and their p-values

textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(",
                   signif(moduleTraitPvalue1, 1), ")", sep = "");
dim(textMatrix1) = dim(moduleTraitCor1)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor1,xLabels = names(datTraits1),
               yLabels = names(MEs1),ySymbols = names(MEs1),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix1,
               setStdMargins = FALSE,
               cex.text = 0.5,zlim = c(-1,1),
               main = paste("Module-Development stage relationships"))

#Gene relationship to trait and important modules: Gene Significance and ModuleMembership
#We quantify associations of individual genes with our trait of interest  by defining Gene Significance GS as(the absolute value of) 
#the correlation between the gene and the trait.  For each module, we also define a quantitative measure of module membership MM
#as the correlation of the module eigengene and the gene expression profile.  This allows us to quantify the similarity of all genes on the array 
#to every module


#for procyclic promastigote ... for other two development stage make the same 
#define the variable of interest
pro= as.data.frame(datTraits1$PRO)
names(pro) = "procyclic"
meta= as.data.frame(datTraits1$META)
names(meta) = "metacyclic"
ama= as.data.frame(datTraits1$AMA)
names(ama) = "amastigote"




# names (colors) of the modules
modNames1 = substring(names(MEs1), 3)
geneModuleMembership1 = as.data.frame(cor(datExpr1, MEs1, use = "p"))
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples1))
names(geneModuleMembership1) = paste("MM", modNames1, sep="");
names(MMPvalue1) = paste("p.MM", modNames1, sep="");

#procyclic
geneTraitSignificance1 = as.data.frame(cor(datExpr1, pro, use = "p"));
GSPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance1), nSamples1));
names(geneTraitSignificance1) = paste("GS.", names(pro), sep="");
names(GSPvalue1) = paste("p.GS.", names(pro), sep="");

#METACYCLIC
geneTraitSignificance2 = as.data.frame(cor(datExpr1, meta, use = "p"));
GSPvalue2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance2), nSamples1));
names(geneTraitSignificance2) = paste("GS.", names(meta), sep="");
names(GSPvalue2) = paste("p.GS.", names(meta), sep="");


#amastigote
geneTraitSignificance3 = as.data.frame(cor(datExpr1, ama, use = "p"));
GSPvalue3 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance3), nSamples1));
names(geneTraitSignificance3) = paste("GS.", names(ama), sep="");
names(GSPvalue3) = paste("p.GS.", names(ama), sep="");



# Intramodular analysis:  identifying genes with high GS and MM
#We can identify genes that have a high significance for interest variable as well as high modulemembership in interesting modules

module1 = "yellow"
column1 = match(module1, modNames1);
moduleGenes1 = moduleColors==module1;

module2 = "blue"
column2 = match(module2, modNames1);
moduleGenes2 = moduleColors==module2

module3 = "turquoise"
column3 = match(module3, modNames1);
moduleGenes3 = moduleColors==module3;



sizeGrWindow(30, 30);
par(mfrow = c(1,1));
par(cex = 1)
#par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))

verboseScatterplot(abs(geneModuleMembership1[moduleGenes1, column1]),
                   abs(geneTraitSignificance1[moduleGenes1, 1]),
                   xlab = paste("Module Membership in", module1, "module"),
                   ylab = "GS for Promastigote",main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, 
                   cex.lab = 1.2, 
                   cex.axis = 1.2,
                   col = module1)

sizeGrWindow(30, 30);
par(mfrow = c(1,1));
par(cex = 1)
verboseScatterplot(abs(geneModuleMembership1[moduleGenes2, column2]),
                   abs(geneTraitSignificance2[moduleGenes2, 1]),
                   xlab = paste("Module Membership in", module2, "module"),
                   ylab = "GS for Metacyclic",main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, 
                   cex.lab = 1.2, 
                   cex.axis = 1.2,
                   col = module2)

sizeGrWindow(30, 30);
par(mfrow = c(1,1));
par(cex = 1)
verboseScatterplot(abs(geneModuleMembership1[moduleGenes3, column3]),
                   abs(geneTraitSignificance3[moduleGenes3, 1]),
                   xlab = paste("Module Membership in", module3, "module"),
                   ylab = "GS for Amastigote",main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, 
                   cex.lab = 1.2, 
                   cex.axis = 1.2,
                   col = module3)






#data summary of information 

names(datExpr1)[moduleColors=="grey"]
# Add annotation information and make a new file for each condition of interest

annot1 = read.csv(file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/WGCNA_annotation_files/LbraM2903_annotationfile.csv", sep = '\t');
genes1 = names(datExpr1)
genes2annot = match(genes1, annot1$query)
sum(is.na(genes2annot))

# Create the starting data frame
geneInfo00 = data.frame(geneID = genes1,
                        annotation = annot1$Predicted.name[genes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance1,
                       GSPvalue1)

modOrder1 = order(-abs(cor(MEs1, pro, use = "p")));


# Add module membership information in the chosen order
for (mod1 in 1:ncol(geneModuleMembership1)){
  oldNames1 = names(geneInfo00)
  geneInfo00 = data.frame(geneInfo00, geneModuleMembership1[, modOrder1[mod1]],
                         MMPvalue1[, modOrder1[mod1]])
  names(geneInfo00) = c(oldNames1, paste("MM.", modNames1[modOrder1[mod1]], sep=""),
                       paste("p.MM.", modNames1[modOrder1[mod1]], sep=""))
}
geneOrder1 = order(geneInfo00$moduleColor, -abs(geneInfo00$GS.procyclic))
geneInfo1 = geneInfo00[geneOrder1, ]

write.csv(geneInfo1, file = "~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/WGCNA_annotation_files/LbraM2903_PRO_geneInfo_tritrypVsncRNA.csv")

#VISUALIZING GENE NETWORK
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM1 = 1-TOMsimilarityFromExpr(datExpr1, power = 12);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM1 = dissTOM1^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM1) = NA;

sizeGrWindow(30,30)

TOMplot(plotTOM1, geneTree, moduleColors, main = "Network heatmap plot, all genes")
# Recalculate module eigengenes


#isualizing the network of eigengenes
# Recalculate module eigengenes
MEs2 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
# Isolate  traits
pro2= as.data.frame(datTraits1$PRO)
names(pro2)="procyclic"
ama2 = as.data.frame(datTraits1$AMA);
names(ama2) = "amastigote"
meta2= as.data.frame(datTraits1$META)
names(meta2) = "metacyclic"
# Add trait to existing module eigengenes
MET1 = orderMEs(cbind(MEs2, ama2))
MET2=orderMEs(cbind(MEs2, pro2))
MET3=orderMEs(cbind(MEs2, meta2))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(30,30);
par(cex = 0.9)
plotEigengeneNetworks(MET2, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.9, xLabelsAngle= 90)

# Plot the dendrogram
sizeGrWindow(9,9)
par(mfrow = c(1,2));
cex1 = 1.0
plotEigengeneNetworks(MET3, "Eigengene dendrogram", marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
plotEigengeneNetworks(MET3, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


#export to Cytoscape
# Recalculate topological overlap if needed
TOM1 = TOMsimilarityFromExpr(datExpr1, power = 12);
# Read in the annotation file
#annot1
# Select modules
modules = modNames1[9];

# Select module probes
inModule = is.finite(match(moduleColors, modules));
modProbes = genes1[inModule];
modGenes = annot1$Predicted.name[match(modProbes, annot1$query)];
modTOM1 = TOM1[inModule, inModule];
dimnames(modTOM1) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(modTOM1,
                               edgeFile = paste("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_tritryp/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("~/Desktop/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_tritryp/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
