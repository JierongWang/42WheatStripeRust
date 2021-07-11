############################# Perform gene-wise association test 
##### Create gene-wise genotype matrix

## Read an pre-process input data
tab <- read.table("all_nsVars_sg", sep="\t", header=TRUE, as.is=TRUE, strip.white=TRUE)
geno <- data.frame(tab[,11:(dim(tab)[2])], stringsAsFactors=FALSE)
## Create table with gene location information
gloc <- read.table("Pst78_genomic.gffgene",sep="\t",header=F)
rownames(gloc) <- gloc$V4

## Read PAV (count) table for identification of gene presence/absence polymorphisms
counts <- read.table("all.pav",sep="\t",header=T)
## Re-order count table columns to match genotype table
counts <- counts[, colnames(geno)]

## Extract gene IDs
unid <- unique(tab$GeneID)
tab2 <- c()
## For each gene with variant(s) ...
for (i in 1:length(unid)){
	## Extract gene ID and row indices of all variants linked to this gene
	id <- unid[i]
	ri <- which(tab$GeneID == id)
	## Extract contig ID
	chr <- unique(tab[ri, "CHROM"])
	## Count number of variants in this gene
	cvar <- dim(tab[ri, ])[1]
	## Extract genotype info for all variants in this gene
	gt <- geno[ri, ]
	## Collapse individual variants into one gene-wise genotype
	gtstr <- apply(gt, 2, function(x){paste(x, collapse="_")})
	## Count how often each genotype occurs
	ft <- sort(table(gtstr), decreasing=TRUE)
	## Extract "main" (i.e. most common) genotype
	major <- unlist(strsplit(names(ft)[1], split="_"))
	## Create list of genotype vectors for each isolate
	gtlist <- strsplit(gtstr, split="_")
	## Re-label gene-wise genotypes (according to frequency of occurence)
	t2 <- c()
	for (j in 1:length(gtlist)){
	## Check which variant calls correpond/do not correspond to main genotype or are 'NA's
        mi <- which(gtlist[[j]] == major)
	oi <- which(gtlist[[j]] != major & gtlist[[j]] != "NA")
        ni <- which(gtlist[[j]] == "NA")
	## Re-label gene-wise genotypes according to frequency of occurence
	## Most frequent genotype is labelled "MAIN"
        if(length(mi) > 0.5 * cvar & length(oi) == 0){
        	t2 <- c(t2, "MAIN")
	## Further genotypes are labelled
        } else if(length(oi) >= 1 & length(ni) < 0.5 * cvar){
		## ... "OTHER" for genes with only two genotypes
		if(length(ft) == 2 & gtstr[j] == names(ft)[2]) {
               		t2 <- c(t2, "OTHER")
		## ... "OTHER1/2/..." for genes with more than two genotypes
             	} else if(length(ft) > 2 & gtstr[j] %in% names(ft)[2:length(ft)]){
                	t2 <- c(t2, paste("OTHER", which(names(ft) == gtstr[j])-1, sep=""))
             	}
		## Cases with no/unclear genotype assignment are labelled
         	} else{
             		## ... "MISS" if the gene is truly not existed (based on PAV table)
			if(counts[id, j] == 0) {
				t2 <- c(t2, "MISS")
			## .. NA in all other cases
			}else{
				t2 <- c(t2, NA)
         		}
		}
     	}
     	names(t2) <- colnames(gt)
     	## Create gene entry and add to gene-wise genotype table
	tab2 <- rbind(tab2, c(GeneID=id, CHROM=chr, SNPs=cvar, t2))
}

rownames(tab2) <- NULL
tab2 <- as.data.frame(tab2, stringsAsFactors=FALSE)

## Among genes without called variants ...
cnv.snpDiff <- counts[setdiff(rownames(counts), unid), ]
tab2b <- apply(cnv.snpDiff, 2, function(x){replace(x, as.numeric(x)==0, "MISS")})
tab2b <- apply(tab2b, 2, function(x){replace(x, as.numeric(x)>=1, "NOM")})
## Add contig information for these genes and adjust table format
tab2b <- data.frame(GeneID=rownames(tab2b), CHROM=gloc[rownames(tab2b),"V1"], SNPs=rep("0", dim(tab2b)[1]), tab2b, stringsAsFactors=FALSE)
rownames(tab2b) <- NULL

## Combine genotyping table and presence/absence table
tab2 <- rbind(tab2, tab2b)
tab2 <- tab2[order(tab2$CHROM, tab2$GeneID),]
rownames(tab2) <- NULL

###### Perform per-gene association test
## Extract genotype information
geno2 <- tab2[, 4:(dim(tab2)[2])]
geno3 <- geno2
# tian data need
#colnames(geno3) <- gsub("X","",colnames(geno3))

### Read phenotype matrix
pheno <- read.table("isolates_phenotype.txt", sep="\t", header=TRUE)
pheno$X <- gsub("-", ".", pheno$X)
rownames(pheno) <- pheno$X
pheno1 <- t(pheno)
pheno1 <- pheno1[,colnames(geno3)]
pheno <- as.data.frame(t(pheno1))
pheno <- pheno[,2:dim(pheno)[2]]

## Create list of Yr alleles
plist2 <- vector("list", length(colnames(pheno)))
names(plist2) <- paste(colnames(pheno), "p", sep="_")
## For each variant position ...
for (i in 1:(dim(geno2)[1])){
        ## ... extract genotype information
        tt2 <- as.factor(unlist(geno2[i,]))
        ## For each Yr allele in phenotyping table
        for (j in 1:length(plist2)){
                 ## ... extract phenotype information
                pp <- as.factor(unlist(pheno[,j]))
                ## If there are at least two different valid (not NA) genotypes,
                ## and at least two samples with valid assignments of both genotype and phenotype ...
                if(length(levels(tt2)) >= 2 & length(levels(pp)) >= 2 & length(which(!is.na(pheno[,j]) & !is.na(tt2))) >= 2){
                        ## ... perform Fisher Test to calculate association p value
                        plist2[[j]] <- c(plist2[[j]], fisher.test(tt2, pp)$p.value)
                ## Else ...
                } else {
                        ## ... set p value to NA 
                        plist2[[j]] <- c(plist2[[j]], NA)
                }

        }

}

## Convert list to data.frame
ptab2 <- do.call(cbind.data.frame, plist2)

## Perform FDR (BH) multiple testing correction of p values
patab2 <- apply(ptab2, 2, function(x){p.adjust(x, method="BH")})
colnames(patab2) <- gsub("_p", "_p.adj", colnames(patab2))

## Create and write output table
out_genes <- data.frame(tab2[,1:(dim(tab2)[2])], ptab2, patab2)
rownames(out_genes) <- NULL
write.table(out_genes, "all_FisherGeneWise", sep="\t", row=FALSE, quote=FALSE)

