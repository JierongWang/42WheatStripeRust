#################### Get the genotype before association test
#### Read tables with annotated functional variants, including missense, start loss, stop loss, stop gain, splice, frameshift, conservative and disruptive inframe variants.
var <- read.table("all.snpIndel.functions.vcf", sep="\t", header=TRUE, as.is=TRUE, strip.white=TRUE)
rownames(var) <- paste(var$CHROM, var$POS, sep=":")
#### Extract codon exchange, amino acid exchange and gene IDs 
td <- unlist(lapply(var$INFO, function(x){strsplit(x, "ANN=")[[1]][2]}))
## Gene ID
idd <- unlist(lapply(td, function(x){strsplit(x, "\\|")[[1]][4]}))
##variant type
vat <- unlist(lapply(td, function(x){strsplit(x, "\\|")[[1]][2]}))
## Codon exchange
ned <- unlist(lapply(td, function(x){strsplit(x, "\\|")[[1]][10]}))
## Amino acid exchange
aed <- unlist(lapply(td, function(x){strsplit(x, "\\|")[[1]][11]}))
#### Extract further information of interest from INFO column
## Overall coverage depth (DP)
dpd <- unlist(lapply(var$INFO, function(x){strsplit(x, "DP=")[[1]][2]}))
dpd <- unlist(lapply(dpd, function(x){strsplit(x, ";")[[1]][1]}))
## Frequency of alternate allele (AF)
afd <- unlist(lapply(var$INFO, function(x){strsplit(x, "AF=")[[1]][2]}))
afd <- unlist(lapply(afd, function(x){strsplit(x, ";")[[1]][1]}))
#### Extract genotype info for each isolate and re-format/simplify genotype
GTD <- var[,10:(dim(var)[2])]

gtd <- c()

for (i in 1:(dim(GTD)[2])) {
	gt <- GTD[,i]
	geno <- vector("character", length=length(gt))
	for (j in 1:length(gt)) {
		y <- unlist(lapply(gt[j], function(x){strsplit(x, ":")}))
		## Simplify genotype information for variants
		##  Classify cases without any variant call information as NA
		if (y[1] == "./."){
			geno[j] <- "NA"
		## Indicate all other cases of low sequence coverage as NA
		} else if (as.numeric(y[3]) < 5){
                        geno[j] <- "NA"
		##  Classify cases same as reference call information as REF
		} else if (y[1] != "./."){
			dd <- as.numeric(unlist(strsplit(y[1], "/")))
			if (y[1] == "0/0"){
				geno[j] <- "REF"
			} else if ( dd[1] > 0 & dd[1] == dd[2]) {
				geno[j] <- paste("ALT", dd[1], sep="")
			## Classify heterozygous genotypes as 'het..'
			} else if ( as.numeric(y[3]) >= 5 & dd[1] != dd[2] ) {
				if (dd[1] == 0) {
					geno[j] <- paste("hetR", dd[2], sep="")
				} else if (dd[1] != 0) {
					geno[j] <- paste("het", gsub("/", "", y[1]), sep="")
				} else {
					geno[j] <- "NA"
				}
			}
		}
	
	}
	gtd <- cbind(gtd, geno)
} 
colnames(gtd) <- colnames(GTD)

## Create annotated genotype table for variants
# snpEff
tdd1 <- data.frame(CHROM=var$CHROM, POS=var$POS, REF=var$REF, ALT=var$ALT, QUAL=round(var$QUAL, 2), DP=dpd, AF=afd, CodonChange=ned, VAT=vat, GeneID=idd, gtd, stringsAsFactors=FALSE)
tdd2 <- data.frame(CHROM=var$CHROM, POS=var$POS, REF=var$REF, ALT=var$ALT, QUAL=round(var$QUAL, 2), DP=dpd, AF=afd, CodonChange=ned, VAT=vat, GeneID=idd, gtd, GTD,stringsAsFactors=FALSE)
##### Output table to file
tdd1 <- tdd1[order(tdd1$CHROM, tdd1$POS), ]
tdd2 <- tdd2[order(tdd2$CHROM, tdd2$POS), ]
write.table(tdd1, "all_nsVars_sg", sep="\t", row=FALSE, quote=FALSE)
write.table(tdd2, "all_nsVars_lg", sep="\t", row=FALSE, quote=FALSE)

