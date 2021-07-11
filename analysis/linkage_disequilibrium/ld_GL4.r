####################### Analysis of recombinatin via linkage
#Only need vcfR for reading in original vcf files and doing some filtering
library("vcfR") 
library("poppr") 
library("ggplot2")

#Read in VCF files for GL4 isolates
vcf_GL4 <- read.vcfR("GL4.snpMV.vcf")
#Convert the VCF files into smaller genlight objects that can be used with poppr/adegenet
gl_GL4 <- vcfR2genlight(vcf_GL4)
save(gl_GL4, file="gl_GL4.Rdata")
#Then reload objects for subsequent analyses if need be
load("gl_GL4.Rdata")
#Set ploidy
gl_list <- list(gl_GL4)
lapply(gl_list, ploidy, 2)
#First remove any NAs
toRemoveGL4 <- is.na(glMean(gl_GL4, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemoveGL4) # position of entirely non-typed loci
gl_GL4_rmNA_LD <- gl_GL4[, !toRemoveGL4]
save(gl_GL4_rmNA_LD, file="gl_GL4_rmNA_LD.Rdata")
#Then reload objects for subsequent analyses if need be
load("gl_GL4_rmNA_LD.Rdata")

#First look at the index of association in randomly sampled sites (5,000 SNP sets repeated 100 times)
#This function will return the standardized index of association ("rbarD")
ia_GL4 <- samp.ia(gl_GL4_rmNA_LD, n.snp = 5000L, reps = 100L)
#In simulations using 10 individuals and 169,588 biallelic SNP loci
### No structure (fully admixed pops)
sex <- glSim(10, 169588, ploid=2, LD=T)
### Most-sex
clone_25 <- glSim(10, 169588, n.snp.struc=42397, ploid=2, LD=T)
### Semi-sex 
clone_50 <- glSim(10, 169588, n.snp.struc=84794, ploid=2, LD=T)
### Most-clonal 
clone_75 <- glSim(10, 169588, n.snp.struc=127191, ploid=2, LD=T)
### Structure (clonal pops)
clone_100 <- glSim(10, 169588, n.snp.struc=169588, ploid=2, LD = T)

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
d1 <- data.frame(ia_GL4, rep("GL4_dataset", length(ia_GL4)))
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
frames <- list(as.data.frame(d1), as.data.frame(d2), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5),as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

# A nonparametric multiple comparison of ranks using Kruskal-Wallis was performed to test for differences across mean ranks.
library(agricolae)

# HSD.test
anova.ia <- aov(lm(ia ~ dset, ia.total))
tukey <- HSD.test(anova.ia, "dset", alpha = 0.001)
# Kruskal wallis test
kruskal.test(ia ~ dset, ia.total) #Kruskal-Wallis chi-squared = 582.36, df = 5, p-value < 2.2e-16
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))

#Plot distributions for GL4
#Set theme to be classic for rest of plots
xname <- factor(ia.total$dset, levels = c("GL4_dataset","sexual","clone_25","clone_50","clone_75","clone_100"))

GL4_plot <- ggplot(ia.total, aes(x = xname, y = ia, fill = as.factor(dset))) +
     geom_boxplot() +
     scale_x_discrete(labels = c("GL4", "No","25%","50%", "75%", "100%")) +
     theme(legend.position="none",axis.title = element_text(size=40),axis.text = element_text(size=28,color="black"), axis.ticks=element_line(size=1,color="black"), panel.background=element_rect(fill="white"),panel.grid.major.y =element_line(color="grey",linetype=0),panel.grid.minor.y =element_line(color="grey",linetype=0),panel.grid.minor.x = element_blank(), axis.line = element_line(size=0.5, colour = "black")) +
     labs(x = "", y = expression(italic(bar("r")["d"]))) + scale_fill_manual(values= c("#c0c000","#00c000","#00c0c0","#F08080","#47b8e0","#34314c"))
ggsave(file="GL4.snpMV5000_25boxplot.pdf",plot=byYear_plot,width =8,height=8)


