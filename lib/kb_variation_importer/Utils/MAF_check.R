#!/usr/bin/env Rscript 

args <- commandArgs(trailingOnly = TRUE)
maf_freq <- read.table(args[1], header =TRUE, as.is=T)
png(args[2])
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()


