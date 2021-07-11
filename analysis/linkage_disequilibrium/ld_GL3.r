####################### Analysis of recombinatin via linkage
#Only need vcfR for reading in original vcf files and doing some filtering
library("vcfR")
library("poppr")
library("ggplot2")

#Read in VCF files for GL3 isolates
vcf_GL3 <- read.vcfR("GL3.snpMV.vcf")
#Convert the VCF files into smaller genlight objects that can be used with poppr/adegenet
gl_GL3 <- vcfR2genlight(vcf_GL3)
save(gl_GL3, file="gl_GL3.Rdata")
#Then reload objects for subsequent analyses if need be
load("gl_GL3.Rdata")
#Set ploidy
gl_list <- list(gl_GL3)
lapply(gl_list, ploidy, 2)
#First remove any NAs
toRemoveGL3 <- is.na(glMean(gl_GL3, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemoveGL3) # position of entirely non-typed loci
gl_GL3_rmNA_LD <- gl_GL3[, !toRemoveGL3]
save(gl_GL3_rmNA_LD, file="gl_GL3_rmNA_LD.Rdata")

#Then reload objects for subsequent analyses if need be
load("gl_GL3_rmNA_LD.Rdata") #1,379,292 biallelic SNPs used for analysis

#First look at the index of association in randomly sampled sites (10,000 SNP sets repeated 100 times)
#This function will return the standardized index of association ("rbarD")
ia_GL3 <- samp.ia(gl_GL3_rmNA_LD, n.snp = 5000L, reps = 100L)
#In simulations using 17 individuals and 299,826 biallelic SNP loci
### No structure (fully admixed pops)
sex <- glSim(17, 299826, ploid=2, LD=T)
### Most-sex
clone_25 <- glSim(17, 299826, n.snp.struc=74957, ploid=2, LD=T)
### Semi-sex 
clone_50 <- glSim(17, 299826, n.snp.struc=149913, ploid=2, LD=T)
### Most-clonal 
clone_75 <- glSim(17, 299826, n.snp.struc=224870, ploid=2, LD=T)
### Structure (clonal pops)
clone_100 <- glSim(17, 299826, n.snp.struc=299826, ploid=2, LD = T)

## IA.sex
ia.sex <- samp.ia(sex, n.snp = 5000L, reps = 100L)
## IA.mostsex
ia.clone.25 <- samp.ia(clone_25, n.snp = 5000L, reps = 100L)
## IA.semisex
ia.clone.50 <- samp.ia(clone_50, n.snp = 5000L, reps = 100L)
## IA.mostclone
ia.clone.75 <- samp.ia(clone_75, n.snp = 5000L, reps = 100L)
## IA.clone
ia.clone.100 <- samp.ia(clone_100, n.snp = 5000L, reps = 100L)

# Summarizing data frames
d1 <- data.frame(ia_GL3, rep("GL3_dataset", length(ia_GL3)))
d2 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d3 <- data.frame(ia.clone.50, rep("clone_50", length(ia.clone.50)))
d4 <- data.frame(ia.clone.75, rep("clone_75", length(ia.clone.75)))
d5 <- data.frame(ia.clone.100, rep("clone_100", length(ia.clone.100)))
d6 <- data.frame(ia.clone.25, rep("clone_25", length(ia.clone.25)))
colnames(d1) <- c("ia","dset")
colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")
colnames(d6) <- c("ia","dset")
ia.total <- rbind(d1, d2, d3, d4, d5, d6)

# Normality tests
frames <- list(as.data.frame(d1), as.data.frame(d2), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5), as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

#A nonparametric multiple comparison of ranks using Kruskal-Wallis was performed to test for differences across mean ranks.
library(agricolae)

# HSD.test
anova.ia <- aov(lm(ia ~ dset, ia.total))
tukey <- HSD.test(anova.ia, "dset", alpha = 0.001)
# Kruskal wallis test
kruskal.test(ia ~ dset, ia.total) #Kruskal-Wallis chi-squared = 582.36, df = 5, p-value < 2.2e-16
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))

#Plot distributions for GL3
#Set theme to be classic for rest of plots
xname <- factor(ia.total$dset, levels = c("GL3_dataset","sexual","clone_25","clone_50","clone_75","clone_100"))

GL3_plot <- ggplot(ia.total, aes(x = xname, y = ia, fill = as.factor(dset))) +
     geom_boxplot() +
     scale_x_discrete(labels = c("GL3", "No","25%","50%", "75%", "100%")) +
     theme(legend.position="none",axis.title = element_text(size=40),axis.text = element_text(size=28,color="black"), axis.ticks=element_line(size=1,color="black"), panel.background=element_rect(fill="white"),panel.grid.major.y =element_line(color="grey",linetype=0),panel.grid.minor.y =element_line(color="grey",linetype=0),panel.grid.minor.x = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
     labs(x = "", y = expression(italic(bar("r")["d"]))) + scale_fill_manual(values= c("#c0c000","#00c000","#00c0c0","#F08080","#F08080","#34314c"))
ggsave(file="GL3.snpMV5000_25boxplot.pdf",plot=byYear_plot,width =8,height=8)

