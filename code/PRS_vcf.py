#!/usr/bin/env python
# coding: utf-8

# In[3]:


import argparse
import math
import pandas as pd
import sys
import vcf
import json


GENOTYPE_ACB = "/home/yiw078/Proj_datasets/ACB/all_chr_ACB_corr.tab"
GENOTYPEv_ACB = "/home/yiw078/Proj_datasets/ACB/sort_all_chr_ACB_corr.recode.vcf.gz"

GENOTYPE_CHS = "/home/yiw078/Proj_datasets/CHS/all_chr_CHS_corr.tab"
GENOTYPEv_CHS = "/home/yiw078/Proj_datasets/CHS/sort_all_chr_CHS_corr.recode.vcf.gz"

GENOTYPE_TSI = "/home/yiw078/Proj_datasets/TSI/all_chr_TSI_corr.tab"
GENOTYPEv_TSI = "/home/yiw078/Proj_datasets/TSI/sort_all_chr_TSI_corr.recode.vcf.gz"

GENOTYPE_ACB_o = "/home/yiw078/Proj_datasets/GWAS/"
GENOTYPEv_ACB_o = "/home/yiw078/Proj_datasets/GWAS/sort_chr1_ACB_o.vcf.gz"

CORRBETA = "/home/yiw078/Proj_datasets/GWAS/corrected_beta.txt"

# Load data
#genotype_ACB = pd.read_csv(GENOTYPE_ACB, sep="\t")

snps = pd.read_csv(CORRBETA, sep="\t")

vcfreader_ACB = vcf.Reader(open(GENOTYPEv_ACB, "rb"))
sample_to_score1_ACB = dict(zip(vcfreader_ACB.samples, [0]*len(vcfreader_ACB.samples)))
sample_to_score2_ACB = dict(zip(vcfreader_ACB.samples, [0]*len(vcfreader_ACB.samples)))


# CHS
#genotype_CHS = pd.read_csv(GENOTYPE_CHS, sep="\t")

vcfreader_CHS = vcf.Reader(open(GENOTYPEv_CHS, "rb"))
sample_to_score1_CHS = dict(zip(vcfreader_CHS.samples, [0]*len(vcfreader_CHS.samples)))
sample_to_score2_CHS = dict(zip(vcfreader_CHS.samples, [0]*len(vcfreader_CHS.samples)))


# TSI
#genotype_TSI = pd.read_csv(GENOTYPE_TSI, sep="\t")

vcfreader_TSI = vcf.Reader(open(GENOTYPEv_TSI, "rb"))
sample_to_score1_TSI = dict(zip(vcfreader_TSI.samples, [0]*len(vcfreader_TSI.samples)))
sample_to_score2_TSI = dict(zip(vcfreader_TSI.samples, [0]*len(vcfreader_TSI.samples)))

# ACB_o
vcfreader_ACB_o = vcf.Reader(open(GENOTYPEv_ACB_o, "rb"))
sample_to_score1_ACB_o = dict(zip(vcfreader_ACB_o.samples, [0]*len(vcfreader_ACB_o.samples)))
sample_to_score2_ACB_o = dict(zip(vcfreader_ACB_o.samples, [0]*len(vcfreader_ACB_o.samples)))


# In[4]:


for i in range(snps.shape[0]):
    chrom = snps["#CHROM"].values[i]
    pos = int(snps["POS"].values[i])
    eff_size1 = float(snps["BETA"].values[i])
    eff_size2 = float(snps["corrBETA"].values[i])
    eff_allele = snps["ALT"].values[i]
    rsid = snps["iID"].values[i]
    matches = vcfreader_ACB.fetch("%s:%s"%(chrom, pos))
    for snp in matches: # should only be one
        alleles = list(map(str, snp.alleles))
        for sample in snp:
            num_eff_alleles = sum(map(lambda x: alleles[int(x)]==eff_allele, sample.gt_alleles))
            #print eff_size1*num_eff_alleles, eff_size2*num_eff_alleles
            sample_to_score1_ACB[sample.sample] += eff_size1*num_eff_alleles
            sample_to_score2_ACB[sample.sample] += eff_size2*num_eff_alleles


# In[5]:


for i in range(snps.shape[0]):
    chrom = snps["#CHROM"].values[i]
    pos = int(snps["POS"].values[i])
    eff_size1 = float(snps["BETA"].values[i])
    eff_size2 = float(snps["corrBETA"].values[i])
    eff_allele = snps["ALT"].values[i]
    rsid = snps["iID"].values[i]
    matches = vcfreader_CHS.fetch("%s:%s"%(chrom, pos))
    for snp in matches: # should only be one
        alleles = list(map(str, snp.alleles))
        for sample in snp:
            num_eff_alleles = sum(map(lambda x: alleles[int(x)]==eff_allele, sample.gt_alleles))
            #print eff_size1*num_eff_alleles, eff_size2*num_eff_alleles
            sample_to_score1_CHS[sample.sample] += eff_size1*num_eff_alleles
            sample_to_score2_CHS[sample.sample] += eff_size2*num_eff_alleles


# In[6]:


for i in range(snps.shape[0]):
    chrom = snps["#CHROM"].values[i]
    pos = int(snps["POS"].values[i])
    eff_size1 = float(snps["BETA"].values[i])
    eff_size2 = float(snps["corrBETA"].values[i])
    eff_allele = snps["ALT"].values[i]
    rsid = snps["iID"].values[i]
    matches = vcfreader_TSI.fetch("%s:%s"%(chrom, pos))
    for snp in matches: # should only be one
        alleles = list(map(str, snp.alleles))
        for sample in snp:
            num_eff_alleles = sum(map(lambda x: alleles[int(x)]==eff_allele, sample.gt_alleles))
            #print eff_size1*num_eff_alleles, eff_size2*num_eff_alleles
            sample_to_score1_TSI[sample.sample] += eff_size1*num_eff_alleles
            sample_to_score2_TSI[sample.sample] += eff_size2*num_eff_alleles


# In[1]:


#ACB_o
#for i in range(snps.shape[0]):
    #chrom = snps["#CHROM"].values[i]
    #pos = int(snps["POS"].values[i])
    #eff_size1 = float(snps["BETA"].values[i])
#    eff_size2 = float(snps["corrBETA"].values[i])
#    eff_allele = snps["ALT"].values[i]
#    rsid = snps["iID"].values[i]
#    matches = vcfreader_ACB_o.fetch("%s:%s"%(chrom, pos))
#    for snp in matches: # should only be one
#        alleles = list(map(str, snp.alleles))
#        for sample in snp:
#            num_eff_alleles = sum(map(lambda x: alleles[int(x)]==eff_allele, sample.gt_alleles))
#            #print eff_size1*num_eff_alleles, eff_size2*num_eff_alleles
#            sample_to_score1_ACB_o[sample.sample] += eff_size1*num_eff_alleles
#            sample_to_score2_ACB_o[sample.sample] += eff_size2*num_eff_alleles


# In[ ]:





# In[5]:


bc_gene.head(10)


# In[33]:


len(sample_to_score2)


# In[6]:


vcfreader_ACB.samples


# In[16]:


dict_list = [sample_to_score2_ACB，sample_to_score2_CHS，sample_to_score2_TSI]
with open('demo3.json', mode='w', encoding='utf-8') as f:
    json.dump(dict_list, f)


# In[ ]:




