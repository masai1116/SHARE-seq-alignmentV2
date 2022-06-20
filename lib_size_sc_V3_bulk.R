#!/usr/bin/Rscript

args <- commandArgs()
# print(args)

dir <- args[6]
Name1 <- args[7]
all <- as.integer(args[8])
dedup <- as.integer(args[9])

Df <- data.frame(sample=Name1, total=all, unique=dedup)
Df$dup <- 1 - Df[ ,3]/Df[ ,2]
colnames(Df)[4] <- "duprate"

Df$temp <- 0
colnames(Df)[5] <- "libsize"
# head(Df)


for (i in 1:nrow(Df)){
  if (Df[i,2] != Df[i,3]) {
    Df[i,5] <- round(uniroot(function(x) (x*(1-exp(-Df[i,2]/x))-Df[i,3]),
                               lower = 0, upper = 1000000000, tol = 1e-7)$root)}
  else{
    Df[i,5] <- Df[i,2]
  }
}
head(Df)
# print(paste("lib size of ", Name1, " is :", Df[1,5] ,sep=""))

temp <- paste(Name1,".libsize.csv", sep="")
File <- paste(dir, temp, sep="/")
write.csv(Df, File, quote=F)

# C/X = 1 - exp( -N/X ) 
# 
# where X = number of distinct molecules in library 
# N = number of read pairs
# C = number of distinct fragments observed in read pairs"

