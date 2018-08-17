require("ggplot2")

args <- commandArgs(trailingOnly = TRUE)

hwe_data<-read.table (file=args[1], header=TRUE)
df <- data.frame("HWE" = hwe_data[,9])
png(args[2])
hwe <- ggplot(df, aes(x=HWE)) + geom_histogram(color="midnightblue", fill="goldenrod2")
print(hwe)
dev.off()

hwe_zoom<-read.table (file=args[3], header=TRUE)
df2 <- data.frame("HWE_ZOOM" = hwe_zoom[,9])
png(args[4])
hwe2 <- ggplot(df2, aes(x=HWE_ZOOM)) + geom_histogram(color="blue", fill="goldenrod")
print(hwe2)
dev.off()
