#!/usr/bin/Rscript
# plot umi/cell or gene/cell

args <- commandArgs(); # print(args)
dir <- args[6]
Name <- args[7]

suppressMessages(library(matrixStats))
suppressMessages(library(reshape2))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))

# load and convert to count matrix
outFile1 <- paste(Name, ".UMIcounts.csv.gz", sep="")
inFile1 <- paste(Name, ".cutoff.bed.gz", sep="")

print("Loading linear gene table...")
linear <- fread(paste(dir, inFile1, sep="/"), header = F, sep="\t")
print("Finished loading file")


if (length(unique(linear$V1)) < 20000){
   print("Grouping UMIs...")
   df <- linear %>% group_by(V1, V2) %>% select(V3) %>% summarise_all(sum)
   df

   print("Converting matrix from long format to wide format...")
   df2 <- pivot_wider(df, names_from = V1, values_from = V3, values_fill = list(V3 = 0))
   df2[1:5, 1:5]

   Df2 <- as(as.matrix(df2[ ,2:ncol(df2)]), "sparseMatrix")
   rownames(Df2) <- df2$V2
   # Df2
} else {
  cells <- unique(linear$V1)
  genes <- unique(linear$V2)
  step=5000

  print(paste("No. genes x cells = ", length(genes), "x", length(cells), sep=""))

  if (length(cells) > step){
      window <- data.frame(start=seq(1, length(cells), step), end=c(seq(step, length(cells), step), length(cells)))
      for (i in 1:nrow(window)){
      	  print(window[i, ])
      	  sub <- linear %>% filter(V1 %in% cells[window[i,1]:window[i,2]])
      	  Mx.sub <- as(acast(sub, V2~V1, value.var="V3", fill=0, fun.aggregate = mean), "sparseMatrix")
      	  emptygenes <- genes[!genes %in% rownames(Mx.sub)]
      	  emptyMx <- as(matrix(data = 0, ncol = ncol(Mx.sub), nrow = length(emptygenes)), "sparseMatrix")
      	  colnames(emptyMx) <- colnames(Mx.sub)
      	  rownames(emptyMx) <- emptygenes
      	  Mx.sub <- rbind(Mx.sub,emptyMx)
      	  if (i ==1){
             Df2 <- Mx.sub
      	   } else{
           Df2 <- cbind(Df2, Mx.sub[rownames(Df2), ])
      	   }
   	}
     } else {
     Df2 <- as(acast(linear, V2~V1, value.var="V3", fill=0, fun.aggregate = mean), "sparseMatrix")
  }
}
print(paste("No. genes x cells = ", nrow(Df2), "x", ncol(Df2), sep=""))

# print(paste("Output to ", outFile1, sep=""))
# write.table(Df2, gzfile(paste(dir, outFile1, sep="/")), sep = "\t", col.names=T, row.names=T, quote=F)

## make h5
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("rhdf5"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("Matrix"))
rna.count <- Df2
print(paste("Writing to ", dir, "/", Name, ".gene.bc.matrices.h5", sep=""))


write10xCounts(
  path = paste(dir, "/", Name, ".gene.bc.matrices.h5", sep=""),
  x = rna.count,
  barcodes = colnames(rna.count),
  gene.id = rownames(rna.count),
  gene.symbol = rownames(rna.count),
  overwrite = FALSE
)

print(paste("Finished making h5 file for ", Name, sep=""))

## make plots      
if (ncol(Df2) <= 50000){
   Df2 <- as.matrix(Df2)
   umiSum <- colSums(Df2); 
   geneSum <- as.data.frame(colCounts(Df2>0))
   rownames(geneSum) <- colnames(Df2)

   print("top10 UMI counts: ")
   head(sort(umiSum, decreasing = T), 10)

   # set cutoff of gene
   Cutoff <- 0
   Idx <- geneSum > Cutoff; # sum(Idx) # number of cell pass filter
   # mean(geneSum[Idx])
   colnames(geneSum) <- "Genes"
   Df3 <- as.data.frame(sort(geneSum$Genes, decreasing = T))
   Df3$Count <- c(1: nrow(Df3)); colnames(Df3) <- c("Gene", "Count")

   file3 <- paste(Name,'.detected.genes.pdf', sep="")
   pdf(paste(dir,file3, sep="/"))
plot(Df3$Count,log10(Df3$Gene), xlab="Barcode rank", ylab = "log10 (Genes)", main = "Detected Genes per Cell", col="darkblue", pch=16)
garbage <- dev.off()

# set cutoff of UMI
Cutoff2 <- 0
Idx2 <- umiSum > Cutoff2; # sum(Idx) # number of cell pass filter 
# mean(umiSum[Idx2])
Df4 <- as.data.frame(sort(umiSum, decreasing = T))
Df4$Count <- c(1: nrow(Df4)); colnames(Df4) <- c("Umi", "Count")

## UMIs per cell
# head(umiSum); length(umiSum)
umiSum <- as.data.frame(umiSum)
geneSum <- as.data.frame(geneSum)

# plot umi vs genes
Df <- cbind(geneSum$Genes, umiSum[match(rownames(geneSum), rownames(umiSum)), ]); colnames(Df) <- c("Genes","UMIs")
# head(Df)
Df <- as.data.frame(Df)

file7 <- paste(Name,'.UMIvsGenes.pdf', sep="")
pdf(paste(dir,file7, sep="/"))
plot(Df$UMIs,Df$Genes, xlab="UMIs", ylab = "Detected Genes", main = "Detected Genes per Cell", col="darkblue", pch=16)
legend("topleft", c(paste("Number of Cells that detected > 0 genes:       ",sum(Df$Genes > 0),sep = ""), 
paste("Number of Cells that detected > 10 genes:     ",sum(Df$Genes > 10), sep = ""),
paste("Number of Cells that detected > 100 genes:   ",sum(Df$Genes > 100),sep = ""),
paste("Number of Cells that detected > 500 genes:   ",sum(Df$Genes > 500),sep = ""),
paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = "")), bty="n") 
garbage <- dev.off()

print(paste("Number of Cells that detected > 0 genes:    ",sum(Df$Genes > 0), sep = ""))
print(paste("Number of Cells that detected > 10 genes:   ",sum(Df$Genes > 10), sep = ""))
print(paste("Number of Cells that detected > 100 genes:  ",sum(Df$Genes > 100), sep = ""))
print(paste("Number of Cells that detected > 500 genes:  ",sum(Df$Genes > 500), sep = ""))
print(paste("Number of Cells that detected > 1000 genes: ",sum(Df$Genes > 1000), sep = ""))
} else {
  print("More than 50k barcodes detected, skip ploting Gene_per_cell...")
}
