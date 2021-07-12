############################ Aligment reads of each sample to reforence genome
# bwa (version 0.7.13-r1126)
# samtools (version 1.8)
# qualimap (version 2.2)
## Creat index of refrence genome
bwa index Pst78_genomic.fna -p Pst78_genomic.fna
## align
bwa aln -t 4 Pst78_genomic.fna GS-15_1_clean.fq.gz > GS-15_1.sai
bwa aln -t 4 Pst78_genomic.fna GS-15_2_clean.fq.gz > GS-15_2.sai
bwa sampe -r '@RG\tID:GS-15\tLB:GS-15\tPL:ILLUMINA\tSM:GS-15' Pst78_genomic.fna GS-15_1.sai GS-15_2.sai GS-15_1_clean.fq.gz GS-15_2_clean.fq.gz | samtools view -Sbh - > GS-15.aln.bam
## sort
samtools sort -@ 4 GS-15.aln.bam -o GS-15.sort.aln.bam
## remove duplicated reads
samtools rmdup GS-15.sort.aln.bam GS-15.sort.rmdup.aln.bam
samtools index -@ 4 GS-15.sort.rmdup.aln.bam
## stat mapping information
qualimap bamqc -bam GS-15.sort.rmdup.aln.bam -outdir GS-15.sort.rmdup.aln.bamqc_result -outformat HTML --java-mem-size=30G
############################# call variants using gatk
# GATK (version 4.0)
## Creat index of refrence genome
java -Xmx30g -jar picard.jar CreateSequenceDictionary R=Pst78_genomic.fna O=Pst78_genomic.fna.dict
samtools faidx Pst78_genomic.fna
## perform gatk
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar HaplotypeCaller -R Pst78_genomic.fna -I GS-15.sort.rmdup.aln.bam -ERC GVCF -O GS-15.gatk.g.vcf.gz
## combine gvcf files of all samples 
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar CombineGVCFs -R Pst78_genomic.fna --variant Tur-1.gatk.g.vcf.gz --variant XZ-21.gatk.g.vcf.gz --variant XZ-20.gatk.g.vcf.gz --variant Chi-1.gatk.g.vcf.gz --variant Ger-1.gatk.g.vcf.gz --variant Eth-1.gatk.g.vcf.gz --variant Eth-2.gatk.g.vcf.gz --variant GS-1.gatk.g.vcf.gz --variant GS-12.gatk.g.vcf.gz --variant GS-14.gatk.g.vcf.gz --variant GS-15.gatk.g.vcf.gz --variant GS-17.gatk.g.vcf.gz --variant GS-3.gatk.g.vcf.gz --variant GS-4.gatk.g.vcf.gz --variant GS-5.gatk.g.vcf.gz --variant GS-8.gatk.g.vcf.gz --variant GS-9.gatk.g.vcf.gz --variant GZ-1.gatk.g.vcf.gz --variant Geo-1.gatk.g.vcf.gz --variant HB-1.gatk.g.vcf.gz --variant Mex-2.gatk.g.vcf.gz --variant GS-20.gatk.g.vcf.gz --variant Pak-1.gatk.g.vcf.gz --variant QH-1.gatk.g.vcf.gz --variant QH-2.gatk.g.vcf.gz --variant QH-3.gatk.g.vcf.gz --variant QH-4.gatk.g.vcf.gz --variant SC-1.gatk.g.vcf.gz --variant SC-2.gatk.g.vcf.gz --variant Tur-2.gatk.g.vcf.gz --variant US-1.gatk.g.vcf.gz --variant XZ-22.gatk.g.vcf.gz --variant XJ-1.gatk.g.vcf.gz --variant XZ-1.gatk.g.vcf.gz --variant XZ-11.gatk.g.vcf.gz --variant XZ-13.gatk.g.vcf.gz --variant XZ-14.gatk.g.vcf.gz --variant XZ-2.gatk.g.vcf.gz --variant XZ-3.gatk.g.vcf.gz --variant XZ-5.gatk.g.vcf.gz --variant XZ-7.gatk.g.vcf.gz --variant YN-2.gatk.g.vcf.gz --output all.gatk.combined.g.vcf.gz
## get vcf
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar GenotypeGVCFs -R Pst78_genomic.fna -O all.gatk.combined.raw.vcf -V all.gatk.combined.g.vcf.gz
## call SNP
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar SelectVariants -R Pst78_genomic.fna -V all.gatk.combined.raw.vcf --select-type-to-include SNP -O all.gatk.snp.raw.vcf
## filter SNP
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar VariantFiltration -R Pst78_genomic.fna -V all.gatk.snp.raw.vcf -O all.gatk.snp.filterM.vcf --filter-expression "FS > 10.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || SOR > 3.0 || MQRankSum < -12.5 || QD < 2.0" --filter-name "MsnpFilter" --cluster-size 3 --cluster-window-size 10
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar SelectVariants -R Pst78_genomic.fna -V all.gatk.snp.filterM.vcf --exclude-filtered -O all.gatk.snp.vcf
vcftools --vcf all.gatk.snp.vcf --max-alleles 2 --min-alleles 2 --min-meanDP 3 --max-meanDP 120 --recode --recode-INFO-all -c > all.snp.vcf
## call indel
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar SelectVariants -R Pst78_genomic.fna -V all.gatk.combined.raw.vcf --select-type-to-include INDEL -O all.gatk.indel.raw.vcf
## filter indel
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar VariantFiltration -R Pst78_genomic.fna -V all.gatk.indel.raw.vcf -O all.gatk.indel.filterM.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "MsnpFilter"
java -Xmx50 -jar gatk-package-4.0.0.0-local.jar SelectVariants -R Pst78_genomic.fna -V all.gatk.indel.filterM.vcf --exclude-filtered -O all.gatk.indel.vcf
vcftools --vcf all.gatk.indel.vcf --max-alleles 2 --min-alleles 2 --min-meanDP 3 --max-meanDP 120 --recode --recode-INFO-all -c > all.indel.vcf
# only retain 1-50bp indel
perl filterIndel.pl all.indel.vcf > all.indel.vcf1
mv all.indel.vcf1 all.indel.vcf
################################### absence/presence variation (PAV)
bedtools genomecov -ibam GS-15.sort.rmdup.aln.bam -bga > GS-15.bedgraph
bedtools intersect -wb -b Pst78_genomic.gffgene.site -a GS-15.bedgraph > GS-15.bedgraph.overlap
perl statPAV.pl GS-15.bedgraph.overlap > GS-15.pav
