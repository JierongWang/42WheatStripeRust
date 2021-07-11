library(vcfR)
library(adegenet)
vcf1 <- read.vcfR("all.snpMV.vcf")
gl1 <- vcfR2genlight(vcf1)
######
grp <- find.clusters(gl1, max.n=10, n.pca=600, scale=FALSE)
##
grp <- find.clusters(gl1, max.n=10, n.pca=600, scale=FALSE, n.clust=4)
grp4 <- as.data.frame(grp$grp)
bic <- as.data.frame(grp$Kstat)
write.table(bic,"bic",sep="\t")
write.table(grp4,"all.snpMV_DAPC_grp4.txt", sep="\t",row=TRUE,col=FALSE,quote=FALSE)
######## plot for BIC
bic <- read.table("bic",header=T)
pic <- ggplot(data=bic,aes(x=K,y=BIC)) + geom_point(color="blue",size=4) +  geom_line(color="blue") + labs(x="Number of clusters",y="BIC") + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) + theme(axis.title = element_text(size=28), axis.text = element_text(size=24,color="black"), axis.ticks=element_line(size=1,color="black"), panel.background=element_rect(fill="white"),panel.grid.major.y =element_line(color="grey",linetype=0),panel.grid.minor.y =element_line(color="grey",linetype=0),panel.grid.minor.x = element_blank(), axis.line = element_line(size=0.5, colour = "black"))
ggsave(file="all_BIC.pdf",plot=pic,width =5,height=4.8)
#####
gl1$pop <- factor(grp$grp)
dapc1 <- dapc(gl1, n.pca=600, n.da=15)
##
scatter(dapc1, grp=dapc1$grp,scree.da=TRUE,posi.da="bottomright",scree.pca=FALSE,legend=TRUE,label.inds=TRUE,bg="white",pch=17:22,cstar=0,col=c("#00b050","#9400d3","#0000c8","#ff8000"))
dapc1_tab <- dapc1$tab
write.table(dapc1_tab,"all.snpMV_DAPC_PCA.txt", sep="\t",row=FALSE,col=TRUE,quote=FALSE)
pca_eig <- as.data.frame(dapc1$pca.eig)
write.table(pca_eig,"all.snpMV_DAPC_pca_eig.txt", sep="\t",row=FALSE,col=TRUE,quote=FALSE)
###################### calculate pairwise Prevosti’s distances
library(poppr)
gn1 <- vcfR2genind(vcf1)
gn1dis <- diss.dist(gn1)
gn1dis_prevosti <- prevosti.dist(gn1)
