library(pheatmap)
################################# pheatmap Av26/V26
library(pheatmap)
vir42 <- read.table("42pheno01_heatmap.txt",header=T)
rownames(vir42) <- vir42$New
hap1 <- vir42[,3:dim(vir42)[2]]
######
annotation_row = data.frame(Av26V26 = factor(rep(c("Av26","V26"), c(31,11))))
rownames(annotation_row) <- rownames(hap1)
ann_colors = list(Av26V26=c(Av26="#ff7473",V26="#47b8e0"))
pheatmap(hap1,cluster_rows=TRUE,cluster_cols=FALSE,show_rownames=TRUE,show_colnames=TRUE,angle_col ="90",fontsize_col=12,treeheight_row=50,legend=FALSE,annotation_row=annotation_row,annotation_legend=TRUE,annotation_colors=ann_colors,color=c("#00C957","yellow"))
######
42N_virPhe.pdf
(9,8)
################################# violin
library(ggplot2)
library(ggpubr)
###
v26 <- read.table("V26.txt",header=T)
mycompar <- list(c("V26","Av26"))
##
mplot <- ggplot(data=v26, aes(x=V26,y=Pheno1,col=V26)) + geom_violin() + geom_boxplot(width=0.1) + labs(x="",y="") + scale_y_continuous(breaks=c(4,6,8,10,12,14,16)) + theme(axis.title = element_text(size=20), axis.ticks=element_line(size=1,color="black"),axis.text.y=element_text(size=22,color="black"), axis.text.x =element_text(size=22,color="black"), panel.background=element_rect(fill="white"),panel.grid.major.y =element_line(color="grey",linetype=0),panel.grid.minor.y =element_line(color="grey",linetype=0),panel.grid.minor.x = element_blank(),legend.position ="none",legend.title=element_blank(),legend.background = element_blank(),legend.key = element_blank(), legend.text=element_text(color="black",size=14), axis.line = element_line(size=0.5, colour = "black")) + stat_compare_means(comparisons = mycompar,method ="wilcox.test",label = "p.format", size=6) + scale_color_manual(values= c("#ff7473","#47b8e0"))
ggsave(file="42N_v26compare.pdf",plot=mplot,width =5,height=5)
