library(dplyr)
library(smfishHmrf)
print("Loaded smfishHmrf")
mem_file <- commandArgs(trailingOnly = TRUE)[1]
nei_file <- commandArgs(trailingOnly = TRUE)[2]
block_file <- commandArgs(trailingOnly = TRUE)[3]
centroid_file <- commandArgs(trailingOnly = TRUE)[4]
beta <- commandArgs(trailingOnly = TRUE)[5]
beta_increment <- commandArgs(trailingOnly = TRUE)[6]
beta_num_iter <- commandArgs(trailingOnly = TRUE)[7]
par_k <- commandArgs(trailingOnly = TRUE)[8]
outdir <- commandArgs(trailingOnly = TRUE)[9]
prefix <- commandArgs(trailingOnly = TRUE)[10]
tolerance <- commandArgs(trailingOnly = TRUE)[11]

beta<-as.double(beta)
beta_increment<-as.double(beta_increment)
beta_num_iter<-as.integer(beta_num_iter)
par_k <- as.integer(par_k)
tolerance <- as.double(tolerance)

#file reading
y<-read.table(mem_file, header=F, row.names=1)
y<-as.matrix(y)
nei<-read.table(nei_file, header=F, row.names=1)
colnames(nei)<-NULL
rownames(nei)<-NULL
nei<-as.matrix(nei)
blocks<-read.table(block_file, header=F, row.names=1)
blocks<-c(t(blocks))
maxblock <- max(blocks)
blocks<-lapply(1:maxblock, function(x) which(blocks == x))
numnei<-apply(nei, 1, function(x) sum(x!=-1))
centroid<-read.table(centroid_file, header=F, row.names=1)
centroid<-as.matrix(centroid)
#parameter setting
k<-par_k
m<-dim(y)[2]
sigma <-array(0, c(m,m,k))
for(i in 1:k){
	sigma[, ,i] <- cov(y)
	print(rcond(sigma[,,i]))
}
mu<-array(0, c(m,k))
kk2<-centroid
for(i in 1:k){
	mu[,i] <- kk2[i,]
}
numcell<-dim(y)[1]
kk_dist<-array(0, c(numcell, k))
for(i in 1:numcell){
	for(j in 1:k){
		kk_dist[i,j] <- dist(rbind(y[i,], mu[,j]), method="euclidean")
	}
}
clust_mem<-apply(kk_dist, 1, function(x) which(x==min(x))) 
lclust<-lapply(1:k, function(x) which(clust_mem == x))
damp<-array(0, c(k))
for(i in 1:k){
	sigma[, , i] <- cov(y[lclust[[i]], ])
	#default tolerance is 1e-60
	di<-findDampFactor(sigma[,,i], factor=1.05, d_cutoff=tolerance, startValue=0.0001)
    if(is.null(di)){
        damp[i] = 0
    }else{
        damp[i] = di
    }
}

#needs y, nei, beta, numnei, blocks, mu, sigma, damp
do_one <- function(prefix, outdir, #all strings
                   par_k, par_y, par_nei, par_beta, par_numnei, par_blocks, 
                   par_mu, par_sigma, par_damp){ #beta is double, par_k is integer

out_file <- sprintf("%s/%s.%.1f.prob.txt", outdir, prefix, par_beta) #hmrfem probability
out_file_2 <- sprintf("%s/%s.%.1f.centroid.txt", outdir, prefix, par_beta) #hmrfem centroids
out_file_3 <- sprintf("%s/%s.%.1f.hmrf.covariance.txt", outdir, prefix, par_beta) #hmrfem covariance
out_file_unnorm <- gsub("prob", "unnormprob", out_file)

tc.hmrfem<-smfishHmrf.hmrfem.multi(y=par_y, neighbors=par_nei, beta=par_beta, numnei=par_numnei, 
blocks=par_blocks, mu=par_mu, sigma=par_sigma, verbose=T, 
err=1e-7, maxit=50, dampFactor=par_damp)

write.table(tc.hmrfem$prob, file=out_file, sep=" ", quote=F, col.names=F, row.names=T)
write.table(tc.hmrfem$unnormprob, file=out_file_unnorm, sep=" ", quote=F, col.names=F, row.names=T)
write.table(t(tc.hmrfem$mu), file=out_file_2, sep=" ", quote=F, col.names=F, row.names=T)
write.table(tc.hmrfem$sigma[,,1], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T)
for(i in 2:par_k){
	write.table(tc.hmrfem$sigma[,,i], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T, append=T)
}
}

beta_current <- beta
for(bx in 1:beta_num_iter){
	print(sprintf("%.3f", beta_current))
	do_one(prefix, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)
	beta_current <- beta_current + beta_increment
}
	