#!/usr/bin/env Rscript 

require("ggplot2")

args <- commandArgs(trailingOnly = TRUE)
maf_freq <- read.table(args[1], header =TRUE, as.is=T)
df <- data.frame("MAF" = maf_freq[, 5])
png(args[2])
maf <- ggplot(df, aes(x=MAF)) + geom_histogram(color="midnightblue", fill="goldenrod2")

print(maf)
dev.off()
# ggsave(
#     args[2],
#     device = "png",
#     width = 3.25,
#     height = 3.25,
#     dpi = 1200
# )
# ggsave(args[2])

