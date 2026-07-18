# running factes for the first time on the sample from 001-101
library(facets)
# testing args
# inf <- '~/local_data/data_store/ctp_102_calling/output/04_FACETS/001-103_Dx.csv.gz'
# outf <- '~/facets_test.pdf'

inf <- snakemake@input[['csv']]
outf <- snakemake@output[['pdf']]
outf.txt <- snakemake@output[['txt']]
normal.samples <- length(snakemake@params[['normal_samples']])
tumor.samples <-length(snakemake@params[['tumor_samples']])

print(paste('NORMAL SAMPLES: ', snakemake@params[['normal_samples']]))
print(paste('NORMAL SAMPLES: ', normal.samples))
print(paste('TUMOR SAMPLES: ', snakemake@params[['tumor_samples']]))
print(paste('TUMOR SAMPLES: ', tumor.samples))

# rcmat <- readSnpMatrix(inf,)
# print(paste('dim(rcmat):', dim(rcmat)))
err.thresh <- Inf
del.thresh <- Inf
pileup <- read.csv(inf, stringsAsFactors = FALSE,
                   colClasses = rep(c("character", "numeric", "character",
                                      "numeric"), c(1, 1, 2, 8)))
if (grepl("chr", pileup$Chromosome[1])) {
    pileup$Chromosome <- gsub("chr", "", pileup$Chromosome)
}
ii <- which(pileup$File1E <= err.thresh & pileup$File1D <=
                del.thresh & pileup$File2E <= err.thresh & pileup$File2D <=
                del.thresh)
rcmat <- pileup[ii, 1:2]

for (i in 5:ncol(pileup)){
    pileup[,i] <- as.numeric(pileup[,i])
}

normal.R.index <- ((1:normal.samples)*4)+1
tumor.R.index <- (((normal.samples+1):(tumor.samples+normal.samples))*4)+1
normal.A.index <- ((1:normal.samples)*4)+2
tumor.A.index <- (((normal.samples+1):(tumor.samples+normal.samples))*4)+2

rcmat$NOR.DP <- rowMeans(pileup[ii, normal.R.index, drop=F]) + rowMeans(pileup[ii, normal.A.index, drop=F])
rcmat$NOR.RD <- rowMeans(pileup[ii, normal.R.index, drop=F])

rcmat$TUM.DP <- rowMeans(pileup[ii, tumor.R.index, drop=F]) + rowMeans(pileup[ii, tumor.A.index, drop=F])
rcmat$TUM.RD <- rowMeans(pileup[ii, tumor.R.index, drop=F])


# try attempts with different seeds until something works
seeds <- c(1234, as.integer(rnorm(100,10000,sd = 2000)))
n.attempt <- 1
while(n.attempt < 100){
    print(paste('Attempt:', n.attempt, "seed:", seeds[n.attempt]))
    set.seed(seeds[n.attempt])
    xx <-  preProcSample(rcmat)
    oo <- procSample(xx)
    fit <- tryCatch(emcncf(oo), error = function(e) {NULL})
    if (!is.null(fit)){
        n.attempt <- 101
    } else {
        n.attempt <- n.attempt +1
    }

}

pdf(outf, height=6, width=8)
plotSample(x=oo,emfit=fit)
logRlogORspider(oo$out, oo$dipLogR)
dev.off()

# get tumor purity and ploiody information

# fit$purity
# fit$ploidy
d <- data.frame(purity=fit$purity, ploidy=fit$ploidy)
rownames(d) <- snakemake@params[['patient_tp']]
write.table(d, outf.txt, sep='\t', quote=F, row.names = T, col.names = T)
