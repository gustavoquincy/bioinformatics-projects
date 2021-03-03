library(locfit)
library(lattice)
library(flexmix)
library(mclust)
library(MASS)
library(cluster)
library(tsne)
library(scatterplot3d)
library(gplots)
library(fpc)
library(ranger)
library(dplyr)
library(Seurat)
 
 
Kessler02_data=Read10X(data.dir="/Volumes/NOBISH/PROJECTS/Kessler02-10X/Kessler02/outs/filtered_feature_bc_matrix/Seuratv2.3.4")
 
dim(Kessler02_data)
 
num_cells = dim(Kessler02_data)[2]
 
summary(colSums(Kessler02_data))
# Min. 1st Qu. Median Mean 3rd Qu. Max.
# 506 9090 29198 33078 46494 234046
 
min_num_cells = floor(num_cells^(1/3))
#par(mfrow = c(1,1))
tmp = apply(Kessler02_data, 1, function(x) sum(x&gt;0))
table(tmp&gt;=min_num_cells)
 
keep = tmp&gt;=min_num_cells
tmp=Kessler02_data[keep,]
 
at_least_one = apply(tmp, 2, function(x) sum(x&gt;0))
min_num_genes = min(at_least_one)
 
colnames(x=Kessler02_data) = paste('Kessler02', colnames(x=Kessler02_data), sep=' ')
Kessler02 = CreateSeuratObject(raw.data=Kessler02_data, min.cells=min_num_cells, min.genes=min_num_genes)
 
save(Kessler02, file="Kessler02_PreQC")
 
mito.genes=grep(pattern="^MT-", x=rownames(x=Kessler02@data), value = TRUE)
percent.mito = colSums(Kessler02@raw.data[mito.genes, ])/colSums(Kessler02@raw.data)
Kessler02 = AddMetaData(object=Kessler02, metadata=percent.mito, col.name="percent.mito")
 
pdf(file="vln_QC.pdf")
VlnPlot(object=Kessler02, features.plot=c("nGene", "nUMI","percent.mito"), nCol=3)
dev.off()
 
pdf(file="Geneplot_QC.pdf")
par(mfrow=c(2,2))
GenePlot(object=Kessler02, gene1="nUMI", gene2="percent.mito")
abline(h=0.59, col="red", lty=2, lwd=2)
abline(v=160000, col="red", lty=2, lwd=2)
GenePlot(object=Kessler02, gene="nUMI", gene2="nGene")
abline(v=160000, col="red", lty=2, lwd=2)
GenePlot(object=Kessler02, gene1="nGene", gene2="percent.mito")
abline(h=0.59, col="red", lty=2, lwd=2)
GenePlot(object=Kessler02, gene1="nGene", gene2="nUMI")
abline(h=160000, col="red", lty=2, lwd=2)
par(mfrow=c(1,1))
dev.off()
 
table(Kessler02@meta.data$nUMI&lt;160000 &amp; Kessler02@meta.data$percent.mito&lt;0.59)
 
# 654
 
## These filtering parameters are determined by the Violin Plot created above, and the table above
Kessler02=FilterCells(Kessler02, subset.names = c("nUMI","percent.mito"), low.thresholds=c(-Inf,-Inf), high.thresholds=c(160000,0.59))
 
Kessler02 = NormalizeData(Kessler02)
 
pdf(file="VariableGenes.pdf")
Kessler02 = FindVariableGenes(Kessler02, do.plot=T, x.low.cutoff=0.0125, y.cutoff=1) #maybe set y.cutoff to 1 or 1.5
dev.off()
 
length(x=Kessler02@var.genes)
 
dim(Kessler02@meta.data)
 
Kessler02=ScaleData(Kessler02, vars.to.regress=c("nUMI","percent.mito"))
 
# Notice the switch to object x instead of Kessler02
x=RunPCA(object=Kessler02, pc.genes=Kessler02@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5, pcs.compute=100)
 
VizPCA(object=x, pcs.use = 1:2)
PCAPlot(object=x, dim.1=1, dim.2=2)
PCHeatmap(object=x, pc.use=1, cells.use=500, do.balanced=T, label.columns=F)
PCHeatmap(object=x, pc.use=1:12, cells.use=500, do.balanced=T, label.columns=F, use.full=F)
 
x=JackStraw(object=x, num.pc=100, num.cores=2, do.par=TRUE)
 
pdf(file="JackStrawPlot.pdf")
JackStrawPlot(object=x, PCs=1:40)
dev.off()
 
pdf(file="ElbowPlot.pdf")
PCElbowPlot(object=x, num.pc=100)
dev.off()
 
  
x=RunTSNE(x, reduction.use="pca", dims.use=1:100, do.fast=FALSE, n.iter=999999999, nthreads=8, seed.use=123)
 
pdf("Kessler02_Clusters.pdf")
x=FindClusters(x, reduction.type="pca", dims.use=1:100, resolution=2, save.SNN=T, plot.SNN=TRUE, algorithm=2, n.iter=999999999, n.start=5000, random.seed=123)
dev.off()
 
pdf("Kessler02_tSNE.pdf")
TSNEPlot(x, do.label=T)
dev.off()
 
pdf(file="ClusterTee.pdf")
x=BuildClusterTree(object=x, do.reorder=TRUE, reorder.numeric=T)
dev.off()
 
node.scores = AssessNodes(object=x)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -&gt; node.scores
 
node.scores
# node oobe
#7 20 0.25714286
#10 17 0.22222222
#8 21 0.18518519
#2 13 0.14537445
#5 18 0.14516129
#9 16 0.13432836
#6 19 0.12962963
#3 14 0.08247423
#1 12 0.03149606
#4 15 0.02362205
 
min_oobe = node.scores$oobe[1]
min_oobe_node = node.scores$node[1]
 
x_final = x
 
while(min_oobe &gt; 0.05) {
x_merged = MergeNode(x_final, node.use=min_oobe_node, rebuild.tree=T)
x_final = BuildClusterTree(object=x_merged, do.reorder=T, reorder.numeric=T)
TSNEPlot(object=x_final, do.return=T, do.label=T)
node.scores=AssessNodes(object=x_final)
node.scores[order(node.scores$oobe, decreasing=T),] -&gt; node.scores
node.scores
min_oobe = node.scores$oobe[1]
min_oobe_node = node.scores$node[1]
}
 
pdf(file="Kessler02_FinalTSNE.pdf")
TSNEPlot(x_final, do.label=T)
dev.off()
 
pdf("GFAP_Expression.pdf")
VlnPlot(x_final, features.plot = "GFAP")
dev.off()
 
pdf("GFAP_FeaturePlot.pdf")
FeaturePlot(x_final, features.plot="GFAP", cols.use=c('grey96','blue'), no.legend=F)
dev.off()
 
pdf("post_cluster_qc.pdf")
VlnPlot(object=x_final, features.plot=c("nGene","nUMI","percent.mito"), nCol=3)
dev.off()
 
library(MAST)
x.markers &lt;- FindAllMarkers(object = x_final, do.print = TRUE, print.bar = TRUE, random.seed = 123, only.pos = TRUE, test.use = "MAST")
x.markers %&gt;% group_by(cluster) %&gt;% top_n(2, avg_logFC)
top2 &lt;- x.markers %&gt;% group_by(cluster) %&gt;% top_n(4, avg_logFC)
 
pdf("TopMarkers_FeaturePlot.pdf")
FeaturePlot(x_final, features.plot = top2$gene, cols.use = c("grey96","blue"), nCol = 5)
dev.off()
 
pdf("TopMarkers_ViolinPlot.pdf")
VlnPlot(x_final, features.plot = top2$gene, x.lab.rot = TRUE, nCol = 5)
dev.off()
 
write.table(x.markers, file = "markers.txt", sep = "\t")
 
 
pdf("TopMarkers_Heatmap.pdf")
DoHeatmap(x_final, genes.use=top2$gene, slim.col.label=T, remove.key=T)
dev.off()
 
 
