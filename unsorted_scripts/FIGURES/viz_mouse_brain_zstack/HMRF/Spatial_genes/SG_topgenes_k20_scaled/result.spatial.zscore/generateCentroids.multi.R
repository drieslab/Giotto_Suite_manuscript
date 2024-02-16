par_k <- commandArgs(trailingOnly = TRUE)[1]
par_seed <- commandArgs(trailingOnly = TRUE)[2]
nstart <- commandArgs(trailingOnly = TRUE)[3]
mem_file <- commandArgs(trailingOnly = TRUE)[4]
centroid_file <- commandArgs(trailingOnly = TRUE)[5]
kmeans_file <- commandArgs(trailingOnly = TRUE)[6]

par_k <- as.integer(par_k)
par_seed <- as.integer(par_seed)
nstart <- as.integer(nstart)
if(par_seed!=-1 & par_seed>0){
set.seed(par_seed)
}
y<-read.table(mem_file, header=F, row.names=1)
y<-as.matrix(y)
k<-par_k
m<-dim(y)[2]
kk<-kmeans(y, k, nstart=nstart, iter.max=100)

write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
		