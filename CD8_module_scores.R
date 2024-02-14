library(data.table)
library(Seurat)

setwd("/scratch/groups/ebutcher/menglan/")

## read scRNAseq data (gene expression, meta data, and barcode info) of total immune cells combined from multiple patient samples from multiple datasets
pipe = fread("gene_expression.csv")
meta = fread("meta.csv")
bc = fread("barcodes.csv")


## select T cells only
bc.t = bc$barcodes[which(bc$`B/T/NK/ILC` == "T")]

pipe = data.frame(pipe)
rownames(pipe) = pipe$V1

pipe = pipe[,-1]
pipe.t = t(as.matrix(pipe))[,bc.t]

meta = as.data.frame(meta)
rownames(meta) = meta$barcodes
meta.t = meta[bc.t,]


## create Seurat object with T cells
obj = CreateSeuratObject(pipe.t, meta.data = meta.t)

## split object by original dataset
obj.list = SplitObject(obj, split.by = "dataset")

 
# identify variable features for each dataset independently
obj.list <- lapply(obj.list, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000))

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(obj.list)

anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

combined <- IntegrateData(anchorset = anchors)


DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.1)


## CD8 modules as defined in Corridoni et al Table S2
subsets = list(naive1 = c("LEF1", "CCR7", "SELL", "TCF7"),
              naive2 = c("LEF1", "TCF7", "ZNF683"),
              trm = c("ITGAE", "GZMA", "CD69", "CD160", "ZNF683"),
              memory = c("GPR183", "CCR7", "ZNF683"),   #CCR7 low
              gzmk.effector1 = c("GZMK", "GZMH", "KLRG1"),
              gzmk.effector2 = c("GZMK", "TNFRSF9", "GZMH"),
              fgfbp2.effector = c("FGFBP2", "GZMH", "FCGR3A"),
              mait = c("RORC", "TRAV1.2", "TRAJ33", "ZBTB16", "NCR3"),  #TRAJ33 not in pipeline
              cd8cd4 = c("CD4", "TNFRSF4", "CD8A", "IL22", "IL17A", "IL1R1", "IL1R2", "CCR6"),
              cd8cd4.foxp3 = c("CD4", "FOXP3", "TNFRSF4", "CD8A", "IL1R1", "IL1R2", "CCR6"),
              il26 = c("IL26", "IL23R", "RORC", "KIT", "IL22", "IL17A", "LAYN", "HAVCR2", "NCR3", "TNFRSF9", "CTLA4", "PDCD1", "AHR", "LAG3", "KIR2DL4", "NCR1", "ENTPD1"),
              tyrobp.neg.iel = c("KLRC2", "KLRD1", "KLRC3", "KIR2DL3", "TNFRSF9", "NCR1", "GZMH", "KIR3DL3", "KIR2DL4", "ENTPD1", "ZNF683"),
              tyrobp.pos.iel = c("KLRC2", "KLRD1", "KLRC3", "KIR2DL3", "TYROBP", "NCR1", "CD160", "KIR3DL2", "KIR3DL3", "KIR2DL4", "FCER1G", "ENTPD1", "ZNF683"),
              cell.cycle = c("MKI67", "ZWINT"))


## calculate score for each CD8 module for each cell

DefaultAssay(combined) <- "RNA"
combined = AddModuleScore(combined, features = subsets)


## export umap for visualization of cell subsets
umap = Embeddings(combined, "umap")

colnames(combined@meta.data)[27:40] = names(subsets)


## manually calculate umap just using marker genes
markers = unique(unlist(subsets, use.names=F))

exp.markers = GetAssayData(combined)

exp.markers = exp.markers[which(rownames(exp.markers) %in% markers),]

umap.man = umap(t(as.matrix(exp.markers)))$layout

colnames(umap.man) = c("UMAP_man1", "UMAP_man2")

out = data.frame(combined@meta.data[,26:40], umap, umap.man)

fwrite(out, "Tcell_UMAP_CD8_mod_scores_Corridoni.csv", row.names=T)



