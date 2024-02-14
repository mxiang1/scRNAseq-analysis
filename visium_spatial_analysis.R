library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(reshape2)
library(ggrepel)

setwd("~/Documents/Projects/Spatial transcriptomics data/")


## import gene lists for GPCRs and secreted proteins
gpcr = read.table("../GPCR_gene_list_HGNC.txt", sep="\t", header=T)
gpcr = gpcr$Approved.symbol

secreted = fread("../Human protein atlas/sa_location_Secreted HPA 2793 genes.tsv")
secreted = secreted$Gene


## read 10x visium data & process
lung = Load10X_Spatial("~/Documents/Projects/Spatial transcriptomics data/Human Lung Cancer (FFPE)/", "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5")

plot1 <- VlnPlot(lung, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(lung, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

lung <- SCTransform(lung, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(lung, features = c("CD3E", "CD3D", "CD4", "CD8A", "ITGAE"))
SpatialFeaturePlot(lung, features = c("CCL19", "CCL21", "CCL2", "APLN", "CXCL13", "CXCL12", "CXCL10", "CXCL9", "CXCL3", "CCL20", "CXCL14"))


#### perform linear regression to identify potential interacting GPCR/secreted protein pairs

exp = as.matrix(GetAssayData(lung, "data"))
secreted.gpcr = exp[which(rownames(exp) %in% c(secreted, gpcr)),]


lm.gene = function(gene, type){
  
  ## perform linear regression for all genes encoding GPCRs or secreted proteins against gene of interest
  lm.res = apply(secreted.gpcr[-which(rownames(secreted.gpcr)==gene),], 1, function(x) summary(lm(x~secreted.gpcr[gene,]))$coefficients[2,c(1,4)])
  
  ## select top 100 most significantly correlated ligands if gene of interest if for a GPCR, and vice versa
  if (type=="secreted") {
    partner = t(lm.res[,which(colnames(lm.res) %in% gpcr)])
    top100 = partner[order(partner[,2]),][1:100,]
  } else if (type=="receptor") {
    partner = t(lm.res[,which(colnames(lm.res) %in% secreted)])
    top100 = partner[order(partner[,2]),][1:100,]
  }
  
  up = rownames(top100)[which(top100[,1]>0)]
  down = rownames(top100)[which(top100[,1]<0)]

  list(up = up, down = down, stats = top100)
}
  
lm.apln = lm.gene("APLN", "secreted")


theme_bw = theme(legend.position = "none", panel.border = element_rect(color="black", fill=NA), panel.grid = element_blank(), panel.background = element_blank())


## make plots for genes positively and negatively correlated with gene of interest, respectively
lm.gene.plot = function(gene, lm.out){
  
  gene.plot = t(secreted.gpcr[c(rownames(lm.out$stats), gene),])
  gene.plot = melt(setDT(as.data.frame(gene.plot)), id.vars = gene)
  gene.plot = aggregate(gene.plot$value, list(gene.plot[,1], gene.plot$variable), mean)
  
  gene.plot$direction = "up"
  gene.plot$direction[which(gene.plot$Group.2 %in% lm.out$down)] = "down"
  
  plot.up = ggplot(gene.plot[which(gene.plot$direction=="up"),], aes(x=Group.1, y=x, label=Group.2)) +
    geom_line(aes(color=Group.2)) +
    geom_point(aes(color=Group.2)) +
    ggtitle(paste0(gene, ": positive correlation")) + xlab(gene) + ylab ("mean expression") + 
    theme_bw
  
  plot.down = ggplot(gene.plot[which(gene.plot$direction=="down"),], aes(x=Group.1, y=x, label=Group.2)) +
    geom_line(aes(color=Group.2)) +
    geom_point(aes(color=Group.2)) +
    ggtitle(paste0(gene, ": negative correlation")) + xlab(gene) + ylab ("mean expression") + 
    theme_bw
  
  plot.up + plot.down
}

lm.gene.plot("APLN", lm.apln)


