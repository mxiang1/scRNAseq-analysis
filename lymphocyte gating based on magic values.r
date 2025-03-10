library(data.table)

genes = fread("genes.csv")  ## dataframe with rows as cells and column as genes
meta = fread("meta.csv") ## meta data for genes.csv with same row order


#### gate on T cells ####
meta = data.frame(meta, Tcells = "", Tinnate = "", Tconv = "", CD4.subsets = "", CD8.subsets = "")

## 1) T cells: CD3E>1.5 & PTPRC>1 & CD14<0.1
meta$Tcells[which(genes$CD3E>1.5 & genes$PTPRC>0.1 & genes$CD14<0.1)] = "yes"


## 2) T innate & myeloid: copy labels author_cell_type
meta$Tinnate[which(meta$author_cell_type %in% c("CD16+_NK", "CD56+_NK", "NK"))] = "NK"
meta$Tinnate[which(meta$author_cell_type == "gdT")] = "gdT"
meta$Tinnate[which(meta$author_cell_type == "ILC")] = "ILC"
meta$Tinnate[which(meta$author_cell_type == "MAIT")] = "MAIT"


## 3) Tconv: T cells -> non gdT/ILC/MAIT/NK/NKT -> CD4 CD4>0.3 & CD8A<0.75 & CD8B<0.5 and CD8 CD8A>1.2 & CD4<0.4
meta$Tconv[which(meta$Tcells == "yes" & meta$Tinnate == "" & genes$CD4>0.3 & genes$CD8A<0.75 & genes$CD8B<0.5)] = "CD4"
meta$Tconv[which(meta$Tcells == "yes" & meta$Tinnate == "" & genes$CD8A>1.2 & genes$CD4<0.4)] = "CD8"


## 4) CD4 subsets: Treg Foxp3>0.25 & IL2RA>0.5 -> Treg hi FOXP3>1.5 & IL2RA>1.5 -> Th17 IL17A>0.1 | IL17F>0.1 -> CD4naive CCR7>1.4 & SELL>1.8 & chemRs<0.1 -> Tfh CXCR5>0.2 & TOX>0.1 -> CD4cm CCR7>1.4 & SELL>1.8 -> Th1 TBX21>0.1 & IFNG>0.2-> CD4others
meta$CD4.subsets[which(meta$Tconv == "CD4" & genes$FOXP3>0.25 & genes$IL2RA>0.5)] = "Treg"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets == "Treg" & genes$FOXP3>1.5 & genes$IL2RA>1.5)] = "Treg hi"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets != "Treg" & meta$CD4.subsets != "Treg hi" & genes$FOXP3<0.1 & (genes$IL17A>0.1 | genes$IL17F>0.1))] = "Th17"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets != "Treg" & meta$CD4.subsets != "Treg hi" & meta$CD4.subsets != "Th17" & genes$FOXP3<0.1 & genes$IL2RA<0.1 & genes$CCR7>1.4 & genes$SELL>1.8 & genes$CXCR3<0.1 & genes$CXCR5<0.1 & genes$CXCR6<0.1 & genes$CX3CR1<0.1 & genes$CCR1<0.1 & genes$CCR2<0.1 & genes$CCR3<0.1 & genes$CCR4<0.1 & genes$CCR5<0.1 & genes$CCR6<0.1 & genes$CCR8<0.1 & genes$CCR10<0.1)] = "CD4naive"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets != "Treg" & meta$CD4.subsets != "Treg hi" & meta$CD4.subsets != "Th17" & meta$CD4.subsets != "CD4naive" & genes$FOXP3<0.1 & genes$IL2RA<0.1 & genes$TOX>0.1 & genes$CXCR5>0.2)] = "Tfh"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets != "Treg" & meta$CD4.subsets != "Treg hi" & meta$CD4.subsets != "Th17" & meta$CD4.subsets != "CD4naive" & meta$CD4.subsets != "Tfh" & genes$FOXP3<0.1 & genes$CCR7>1.4 & genes$SELL>1.8)] = "CD4cm"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets != "Treg" & meta$CD4.subsets != "Treg hi" & meta$CD4.subsets != "Th17" & meta$CD4.subsets != "CD4naive" & meta$CD4.subsets != "Tfh" & meta$CD4.subsets != "CD4cm" & genes$FOXP3<0.1 & genes$TBX21>0.1 & genes$IFNG>0.2)] = "Th1"

meta$CD4.subsets[which(meta$Tconv == "CD4" & meta$CD4.subsets == "")] = "CD4 others"


## 5) CD8 subsets: CD8naive CCR7>1.4 & SELL>1.8 & chemRs<0.1 -> CD8cm CCR7>1.4 & SELL>1.8 -> CD8em CCR7<1.4 | SELL<1.8 -> CD8trm CCR7<1.4 | SELL<1.8 & ITGAE>0.8

meta$CD8.subsets[which(meta$Tconv == "CD8" & genes$CCR7>1.4 & genes$SELL>1.8 & genes$CXCR3<0.1 & genes$CXCR5<0.1 & genes$CXCR6<0.1 & genes$CX3CR1<0.1 & genes$CCR1<0.1 & genes$CCR2<0.1 & genes$CCR3<0.1 & genes$CCR4<0.1 & genes$CCR5<0.1 & genes$CCR6<0.1 & genes$CCR8<0.1 & genes$CCR10<0.1)] = "CD8naive"

meta$CD8.subsets[which(meta$Tconv == "CD8" & meta$CD8.subsets != "CD8naive" & genes$CCR7>1.4 & genes$SELL>1.8)] = "CD8cm"

meta$CD8.subsets[which(meta$Tconv == "CD8" & meta$CD8.subsets != "CD8naive" & meta$CD8.subsets != "CD8cm" & genes$ITGAE<0.8)] = "CD8 TEM"

meta$CD8.subsets[which(meta$Tconv == "CD8" & meta$CD8.subsets != "CD8naive" & meta$CD8.subsets != "CD8cm" & genes$ITGAE>0.8)] = "CD8 TRM"


meta$t.subsets = ""
meta$t.subsets[which(meta$Tconv=="CD4")] = meta$CD4.subsets[which(meta$Tconv=="CD4")]
meta$t.subsets[which(meta$Tconv=="CD8")] = meta$CD8.subsets[which(meta$Tconv=="CD8")]


#### gate on B cells ####

meta$b.subsets = ""
meta$b.subsets[which(meta$author_cell_type == "Plasma_B")] = "Plasma"
meta$b.subsets[which(meta$author_cell_type == "atypical_B")] = "B atypical"
meta$b.subsets[which(meta$author_cell_type == "naive_B")] = "B naive"
meta$b.subsets[which(meta$author_cell_type %in% c("IGHMhi_memory_B", "IGHMlo_memory_B"))] = "B memory"
meta$b.subsets[which(meta$b.subsets=="B memory" & genes$FCRL4>1)] = "FCRL4+ MBC"


## select genes to include
goi = c("GPR25", "GPR15", "CCR9", "CCR10", "INAVA")
gating.genes = c("PTPRC", "CD14", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "IL17A", "IL17F", "CCR7", "SELL", "CXCR5", "TOX", "TBX21", "IFNG", "ITGAE")
chem.receptors = c("CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR6", "CX3CR1", "XCL1", "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR8", "CMKLR1")
smoking.genes = c("AHR", "AHRR", "CYP1A1", "CYP1B1", "HIF1A", "ARNT", "EPAS1", "F2RL3")
sg.genes = c("HNRNPLL", "PTPRC", "TCHP", "FLNA", "SPSB2", "CCL4", "IRF5", "TECR", "SH3YL1", "CLEC2D", "POLB", "CD83", "BCL2A1")


genes.include = unique(c(goi, gating.genes, chem.receptors, smoking.genes, sg.genes))

meta = data.frame(meta, genes[, ..genes.include])

fwrite(meta, "AIDA_meta_markers.csv")
