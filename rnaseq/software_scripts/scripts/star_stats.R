files = Sys.glob("02-STAR_alignment/*/*_ReadsPerGene.out.tab")
stats = t(as.data.frame(lapply(files, function(x){read.delim(x,header=F,nrows=4)[,4]})))
samples=read.delim("samples.txt",header=F)
rownames(stats) = samples[,1]
cn = lapply(files, function(x){read.delim(x,header=F,nrows=4)})[[1]][,1]
colnames(stats) = cn
stats = cbind(stats,total_in_feature=rowSums(stats))
stats = cbind(sample=rownames(stats), stats)
write.table(stats, file="summary_star_alignments.txt", col.names=T, row.names=F, quote=F, sep="\t")
