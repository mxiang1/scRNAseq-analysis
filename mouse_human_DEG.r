library(Seurat)
library(openxlsx)

## import human and mouse data

setwd("~/oak/scRNAseq/cellRanger3.0/human/LN4_LEC/outs/")

human = Read10X("filtered_feature_bc_matrix/")
colnames(human) = paste0(colnames(human), "_LN4_LEC")


## use a human to mouse gene conversion list to convert human genes to mouse orthologs
genes = read.csv("gene conversion/mouse_to_human_genes_unique_USETHISFORUMAP.csv")

human = data.frame(human_symbol = rownames(human), human)
human.mousegenes = merge(genes, human, by="human_symbol")

rownames(human.mousegenes) = human.mousegenes$mouse_symbol
human.mousegenes = human.mousegenes[,-c(1:5)]


## create a seurat object for human data
human.obj = CreateSeuratObject(counts = human.mousegenes, min.cells = 3, min.features = 100)
human.obj = NormalizeData(human.obj, normalization.method = "LogNormalize", scale.factor = 10000)



## integrate human and mouse data using Seurat (CCA)

mouse.obj = readRDS("Mouse.rds")

comp = list(mouse.obj, human.obj)

for (i in 1:length(comp)) {
  comp[[i]] <- FindVariableFeatures(comp[[i]], selection.method = "vst",
                                    nfeatures = 2000, verbose = FALSE)
}

anchors <- FindIntegrationAnchors(object.list = comp, dims = 1:30)

comp.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
Idents(comp.integrated) = comp.integrated$subsets


DefaultAssay(comp.integrated) <- "integrated"

comp.integrated <- ScaleData(comp.integrated, verbose = FALSE)
comp.integrated <- RunPCA(comp.integrated, npcs = 30, verbose = FALSE)

comp.integrated = FindNeighbors(comp.integrated, dims = 1:10)
comp.integrated = FindClusters(comp.integrated, resolution = 0.5)

comp.integrated <- RunUMAP(comp.integrated, dims = 1:10)
comp.integrated <- RunTSNE(comp.integrated, dims = 1:10)

DimPlot(comp.integrated, reduction="umap", split.by = "sample")
DimPlot(comp.integrated, reduction="tsne", split.by = "sample", group.by = "subsets")



#### DEG analysis for mouse and human

mouse.allmarkers = FindAllMarkers(mouse.obj, test.use = "negbinom", only.pos = T, logfc.threshold = log(1)) # p<0.01
human.allmarkers = FindAllMarkers(human.obj, test.use = "negbinom", only.pos = T, logfc.threshold = log(1)) # p<0.01


## shared lists between mouse and human for each cluster

shared=list()

for (i in 1:5){
  
  mouse.i = mouse.allmarkers[which(mouse.allmarkers$cluster_mouse==clusters[[i]]),]
  human.i = human.allmarkers[which(human.allmarkers$cluster_human==clusters[[i]]),]
  
  shared[[i]] = merge(mouse.i, human.i, by.x="gene_mouse", by.y="mouse_symbol")
  shared[[i]] = shared[[i]][order(-shared[[i]]$avg_logFC_mouse),]
  shared[[i]] = shared[[i]][!duplicated(shared[[i]]),]
}

names(shared) = c("0 cLEC", "1 fLEC", "3 Ptx3-LEC", "4 Marco-LEC", "7 valve")

write.xlsx(shared, "Human_mouse_shared_DEG_5subsets_FC1.1_p0.01_Nov24.xlsx")


# unique lists - MOUSE

uni=list()

for (i in 1:5){
  
  mouse.i = mouse.allmarkers[which(mouse.allmarkers$cluster_mouse==clusters[[i]]),]
  human.all.i = human.allmarkers[which(human.allmarkers$cluster_human==clusters[[i]]),]
  
  uni[[i]] = merge(mouse.i, human.all.i, by.x="gene_mouse", by.y="mouse_symbol", all.x=T)
}

uni.mouse = uni
names(uni.mouse) = c("0 cLEC", "1 fLEC", "3 Ptx3-LEC", "4 Marco-LEC", "7 valve")

write.xlsx(uni.mouse, "Mouse_unique_DEG_5subsets_FC1.1_p0.01_Nov24.xlsx")


# unique lists - HUMAN

uni=list()

for (i in 1:5){
  
  human.i = human.allmarkers[which(human.allmarkers$cluster_human==clusters[[i]]),]
  mouse.all.i = mouse.allmarkers[which(mouse.allmarkers$cluster_mouse==clusters[[i]]),]
  
  uni[[i]] = merge(human.i, mouse.all.i, by.x="mouse_symbol", by.y="gene_mouse", all.x=T)
}

uni.human = uni
names(uni.human) = c("0 cLEC", "1 fLEC", "3 Ptx3-LEC", "4 Marco-LEC", "7 valve")

write.xlsx(uni.human, "Human_unique_DEG_5subsets_FC1.1_p0.01_Nov24.xlsx")

