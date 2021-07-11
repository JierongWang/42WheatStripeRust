############################ association test per vartiants
##### Read an pre-process input data
tab <- read.table("all_nsVars_sg", sep="\t", header=TRUE, as.is=TRUE, strip.white=TRUE)
### Read genotype data 
geno <- data.frame(tab[,11:(dim(tab)[2])], stringsAsFactors=FALSE)
geno3 <- geno

### Read phenotype matrix
pheno <- read.table("isolates_phenotype.txt", sep="\t", header=TRUE)
pheno$X <- gsub("-", ".", pheno$X)
rownames(pheno) <- pheno$X
pheno1 <- t(pheno)
pheno1 <- pheno1[,colnames(geno3)]
pheno <- as.data.frame(t(pheno1))
pheno <- pheno[,2:dim(pheno)[2]]

##### Perform per-Variant association test 
### Create list of Yr alleles
plist <- vector("list", length(colnames(pheno)))
names(plist) <- paste(colnames(pheno), "p", sep="_")
## For each variant position ...
for (i in 1:(dim(geno)[1])){
        ## ... extract genotype information
        tt <- as.factor(unlist(geno[i,]))
        ## For each Yr allele in phenotyping table
        for (j in 1:length(plist)){
                 ## ... extract phenotype information
                pp <- as.factor(unlist(pheno[,j]))
                ## If there are at least two different valid (not NA) genotypes,
                ## and at least two samples with valid assignments of both genotype and phenotype ...
                if(length(levels(tt)) >= 2 & length(levels(pp)) >= 2 & length(which(!is.na(pheno[,j]) & !is.na(tt))) >= 2){
                        ## ... perform Fisher Test to calculate association p value
                        plist[[j]] <- c(plist[[j]], fisher.test(tt, pp)$p.value)
                ## Else ...
                } else {
                        ## ... set p value to NA 
                        plist[[j]] <- c(plist[[j]], NA)
                }

        }

}

## Convert list to data.frame
ptab <- do.call(cbind.data.frame, plist)


## Perform FDR (BH) multiple testing correction of p values
patab <- apply(ptab, 2, function(x){p.adjust(x, method="BH")})
colnames(patab) <- gsub("_p", "_p.adj", colnames(patab))

## Create and write output table
out_sites <- data.frame(tab[,1:(dim(tab)[2])], ptab, patab)
rownames(out_sites) <- NULL
write.table(out_sites, "all_FisherPerVar", sep="\t", row=FALSE, quote=FALSE)

