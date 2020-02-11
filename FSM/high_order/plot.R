
####################################################################
# plot conductance vs group, cluster size vs group
library(ggplot2)
setwd("/Users/tao/FS2016_10/OurESURemoveD/high_order/Exp/size_5")
data = read.table("results_5_even.size_condc.txt", header=T, sep = "\t", stringsAsFactors=F)

data_10k_50k = subset(data, select=c("X10k_50k.cluster_size","X10k_50k.condc"))
data_1k_5k = subset(data, select=c("X1k_5k.cluster_size","X1k_5k.condc"))
data_100_500 = subset(data, select=c("X100_500.cluster_size","X100_500.condc"))
data_10_50 = subset(data, select=c("X10_50.cluster_size","X10_50.condc"))
colnames(data_10k_50k) =c("ClusterSize", "Condc")
colnames(data_1k_5k) =c("ClusterSize", "Condc")
colnames(data_100_500) =c("ClusterSize", "Condc")
colnames(data_10_50) =c("ClusterSize", "Condc")
data_10k_50k$group = "10k-50k"
data_1k_5k$group = "1k-5k"
data_100_500$group = "100-500"
data_10_50$group = "10-50"

new_data = rbind(data_10k_50k, data_1k_5k, data_100_500)
new_data$group = factor(new_data$group, levels=c("10k-50k","1k-5k","100-500"))
new_data$SCondc = new_data$Condc/new_data$ClusterSize

gg1 = ggplot(new_data, aes(x=group,y=Condc)) + geom_violin() + geom_jitter(width=0.25, na.rm=TRUE) + xlab("Number of subgraph instances") + ylab("Motif conductance")
ggsave("gg1.pdf", width=5,height=5)
gg2 = ggplot(new_data, aes(x=group,y=ClusterSize)) + geom_violin() + geom_jitter(width=0.25, na.rm=TRUE) + xlab("Number of subgraph instances") + ylab("Largest cluster size")
ggsave("gg2.pdf", width=5,height=5)
gg3 = ggplot(new_data, aes(x=group,y=SCondc)) + geom_violin() + geom_jitter(width=0.25, na.rm=TRUE) + xlab("Number of subgraph instances") + ylab("Normalized Motif conductance")
ggsave("gg3.pdf", width=5,height=5)
pdf("NormalizedConduc.pdf", width=5, height=5)
boxplot(SCondc~group, data=new_data, outline=F, ylab="Normalized Conductance")
dev.off()
pdf("ClusterSize.pdf", width=5, height=5)
boxplot(ClusterSize~group, data=new_data, outline=F, ylab="Cluster Size")
dev.off()
####################################################################
# plot cam distribution
data = read.table("cam_distribution.txt", header =F, sep = "\t", stringsAsFactors =F)
pdf("cam_distribution.pdf", width = 5, height = 5)
hist(log10(data$V2),xlab="Number of subgraph instances(log10 transformed)", ylab="Frequency of subgraphs", main="")
dev.off()

####################################################################
# plot cluster coverage
library(pheatmap)

data = read.table("cluster_overlaps.txt", header=T, stringsAsFactors=F, sep = " ", check.names=F)
mat = as.matrix(data)
rownames(mat) = colnames(mat)

pheatmap(mat, cluster_rows=F, cluster_cols=F, fontsize=5, filename = "cluster_overlaps.pdf", width = 7, height = 7)
####################################################################
# plot zscore
data = read.table("/Users/tao/FS2016_10/OurESURemoveD/high_order/Exp/size_5/results_5_even.txt.zscore", header=F, stringsAsFactors=F)
data$group = factor(data$V1, levels=c("10k_50k", "1k_5k","100_500", "10_50"))
pdf("Zscore.pdf",  width=5, height=5)
boxplot(V3~group, data, outline=F, ylab="Z-score")


