setwd("~/Documents/Projects/conservation")

library(seqinr)
library(Peptides)
library(biomaRt)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)


## import all protein sequences for 6 species: human, mouse, rat, rabbit, dog, cow

files = c("Bos_taurus.ARS-UCD1.2.pep.all.fa", "Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa", "Homo_sapiens.GRCh38.pep.all.fa", "Mus_musculus.GRCm39.pep.all.fa", "Oryctolagus_cuniculus.OryCun2.0.pep.all.fa", "Rattus_norvegicus.mRatBN7.2.pep.all.fa")  ## downloaded from ensembl

aa = lapply(files, function(x) read.fasta(paste0("../Ensembl protein seq/", x), seqtype="AA", as.string=T, set.attributes=T))

names(aa) = c("Cow", "Dog", "Human", "Mouse", "Rabbit", "Rat")


## extract raw sequences and annotation information which contains gene name and other information
seq = lapply(aa, function(x) unlist(lapply(x, function(y) y[1])))
id = lapply(aa, function(x) unlist(lapply(x, function(y) attr(y, "Annot"))))

for (i in names(seq)){
  names(seq[[i]]) = id[[i]]
}


## only keep proteins longer than 6 amino acids (since we are comparing the C-terminal 6 AA across species
clean = lapply(seq, function(x) x[nchar(x)>=6])


## format annotation info so that only gene names are kept
clean = lapply(clean, function(x) setNames(x, paste0(gsub(".+gene_symbol:([A-Za-z0-9_./-]+) .+", "\\1", names(x)), gsub(" .+", "", names(x)))))


## extract last 6 AA of each protein and create a data frame with ensembl gene ID, gene name, protein ID, and sequence of the last 6 AA for each species
last6 = lapply(clean, function(x) gsub(".+(.{6})$", "\\1", x))
last6.df = lapply(last6, function(x) unique(data.frame(names(x), unlist(lapply(strsplit(names(x), split = ">"), "[" , 1)), unlist(lapply(strsplit(names(x), split = ">"), "[", 2)), x)))

for (i in names(last6.df)){
  colnames(last6.df[[i]]) = c(paste0(i, ".ID"), paste0(i,".gene.name"), paste0(i, ".Protein.ID"), paste0(i, ".last6AA"))
}


## calculate pI and MW for all proteins and gate on those whose human orthologs have a pI > 8 & MW < 30 kD (we are interested in positively charged smaller chemokines)

human.gate = data.frame(ID = names(clean[["Human"]]), gene = unlist(lapply(strsplit(names(clean[["Human"]]), split = ">"), "[" , 1)), seq = clean[["Human"]], MW = mw(clean[["Human"]]), pI = pI(clean[["Human"]]))
human.gate = human.gate[-grep("*", human.gate$seq, fixed=T),]

human.gate = human.gate[which(human.gate$MW<30000 & human.gate$pI>8),]




#### process gene ortholog table of 6 species exported from Ensembl: human, mouse, rat, rabbit, dog, cow ####

orthologs = read.csv("Ensembl_genes_human_mouse_rat_rabbit_dog_cow.txt")

## remove rows where human gene is empty or only human gene exists
ortho.unique = orthologs[-which(orthologs[,1]=="" | orthologs[,-1]==""),]


## keep only genes whose human orthologs fall in the gate in line 52
ortho.uni.gated = ortho.unique[which(ortho.unique$Gene.name %in% human.gate$gene),]


## merge tables containing the last 6 C-terminal AA from all 6 species
ortho.seq = merge(ortho.uni.gated, last6.df[["Human"]], by.x = "Gene.name", by.y = "Human.gene.name")
ortho.seq = merge(ortho.seq, last6.df[["Mouse"]], by = "Mouse.gene.name", all.x = T)
ortho.seq = merge(ortho.seq, last6.df[["Rat"]], by = "Rat.gene.name", all.x = T)
ortho.seq = merge(ortho.seq, last6.df[["Rabbit"]], by = "Rabbit.gene.name", all.x = T)
ortho.seq = merge(ortho.seq, last6.df[["Dog"]], by = "Dog.gene.name", all.x = T)
ortho.seq = merge(ortho.seq, last6.df[["Cow"]], by = "Cow.gene.name", all.x = T)

saveRDS(ortho.seq, "ortho.seq.rds") # 153,338,581 rows


## import list of secreted proteins
secreted = read.table("../Human protein atlas/sa_location_Secreted HPA 2793 genes.tsv", sep="\t", header=T)

secreted = human.gate[which(human.gate$gene %in% secreted$Gene),]


## keep only secreted proteins with pI > 8 & MW < 30 kD
ortho.seq.se = ortho.seq[which(ortho.seq$Gene.name %in% secreted$gene),]  # 25,799,535 rows

ortho.seq.uni = unique(ortho.seq.se[,c(6:1,9,12,15,18,21,24)])  # 192426 rows

ortho.seq.ids = unique(ortho.seq.se[,c(6:1,9,12,15,18,21,24,8,11,14,17,20,23)]) # 25799535 rows


## save each sequence in a fasta file with gene and species name in the header to be used for BLASTP alignment
fa = data.frame(id = paste(seq(1:nrow(ortho.seq.uni)), ortho.seq.uni$Gene.name, sep="_"), ortho.seq.uni)

species = c("human", "mouse", "rat", "rabbit", "dog", "cow")
colnames(fa)[2:7] = species

fa_header = function(df){
  for (i in 1:nrow(df)){
    for (j in species) {
      if (df[i,j]!="") { df[i,j] = paste0(">", df[i,j], "_", i, "_", j) }
    }
  }
  df
}

fa = fa_header(fa)


## for BLASTP, human sequences are provided as query and sequences of the other species are subject
for (i in 1:nrow(fa)){
  write.table(t(fa[i,c(2,8)]), paste0("fasta/", fa[i,1], ".query.fa"), row.names = F, col.names = F, quote = F, na = "")
  write.table(t(fa[i,c(3,9,4,10,5,11,6,12,7,13)]), paste0("fasta/", fa[i,1], ".subject.fa"), row.names = F, col.names = F, quote = F, na = "")
}


## run and process blastp results in Sherlock (Stanford HPC)

setwd("/scratch/groups/ebutcher/menglan/aa")

hits = read.csv("blastp_summary.csv")

files = hits$file[which(hits$nrow!=5 & hits$n.!=5)]

blastp = lapply(files, function(x) read.table(paste0("blastp/", x)))
names(blastp) = gsub(".out.txt", "", files)
## blastp output columns: "query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"


## sum up all the BLAST scores for each ortholog set
blastp.sum = lapply(blastp, function(x) data.frame(matched=nrow(x), score.sum=sum(x$V12)))
blastp.sum = do.call(rbind, blastp.sum)
blastp.sum = data.frame(id = rownames(blastp.sum), blastp.sum)

blastp.seq = merge(fa, blastp.sum, by = "id", all.x = T)


## calculate the average BLAST score for each ortholog set
n.seq = apply(blastp.seq[,8:13], 1, function(x) 5-sum(is.na(x)))
score.mean = blastp.seq$score.sum/n.seq

blastp.seq = data.frame(gene.name.human = gsub(".+_", "", blastp.seq$id), blastp.seq, n.seq, score.mean)

blastp.seq = blastp.seq[order(-blastp.seq$score.mean),]

fwrite(blastp.seq, "blastp_secretome_MW30kD_pI8_alltx_6species.csv", na = "")



