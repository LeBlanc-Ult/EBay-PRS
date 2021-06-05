#!/bin/sh

snpfile=~/Proj_datasets/rsIDs_new.txt
DATADIR=~/Proj_datasets/${pop}/all_chr_${pop}_merge
OUTDIR=~/Proj_datasets

DATADIR2=/datasets/cs284-sp21-A00-public/1000Genomes
KGVCF=${DATADIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
OUTDIR=~/Proj_datasets

#snp filter
for pop in 'ACB' 'CHS' 'TSI'
do
	for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
	do
        vcftools --gzvcf ${DATADIR2}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --snps ${snpfile} --keep ${OUTDIR}/${pop}.txt --recode --recode-INFO-all --out ${OUTDIR}/${pop}/chr${chrom}_${pop}_subset
    done
done

#merge vcf
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
   	echo $chrom
    vcf-concat -f ${OUTDIR}/${pop}/chrom_1_22_${pop}.txt > ${OUTDIR}/${pop}/all_chr_${pop}_merge
done

#convert to .tab
for pop in 'ACB' 'CHS' 'TSI'
do
    cat ${OUTDIR}/${pop}/all_chr_${pop}_corr.recode.vcf | vcf-to-tab > ${OUTDIR}/${pop}/all_chr_${pop}_corr.tab
done

#sort and index for PRS calculate
for pop in 'ACB' 'CHS' 'TSI'
do
    bgzip ${OUTDIR}/${pop}/all_chr_${pop}_corr.recode.vcf
    bcftools sort ${OUTDIR}/${pop}/all_chr_${pop}_corr.recode.vcf.gz -o ${OUTDIR}/${pop}/sort_all_chr_${pop}_corr.recode.vcf
    bgzip ${OUTDIR}/${pop}/sort_all_chr_${pop}_corr.recode.vcf
    bcftools index -t ${OUTDIR}/${pop}/sort_all_chr_${pop}_corr.recode.vcf.gz
done

#preprocess simulated data for PRS
OFFSPRING="~/Proj_datasets/GWAS/chr1_ACB_offspring.vcf"
OUTDIR=~/Proj_datasets

bgzip ${OFFSPRING}
bcftools sort ${OFFSPRING}.gz -o ${OUTDIR}/GWAS/sort_chr1_ACB_offspring.vcf
bgzip ${OUTDIR}/GWAS/sort_chr1_ACB_offspring.vcf
bcftools index -t ${OUTDIR}/GWAS/sort_chr1_ACB_offspring.vcf.gz



