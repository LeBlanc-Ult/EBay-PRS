#install.packages("sim1000G")
#install.packages("tidyr")
#install.packages("pryr")

require(sim1000G)
require(tidyr)
require(pryr)

readvcf <- function(filename = "data.vcf", thin = NA, maxNumberOfVariants = 400,
                   min_maf = 0, max_maf = NA, region_start = NA, region_end = NA)
  {
  if(grepl(".rdata", filename)) {
    cat("[#.......] Loading VCF environment from rdata file\n")
    
    e1 <- new.env()
    load(filename, envir = e1)
    
    vcf <- e1$vcf
    cat("[##......] Number of variants in VCF: ", nrow(vcf$vcf), "\n");
    cat("[##......] Number of individuals in VCF: ", ncol(vcf$gt1), "\n");
    
    return(vcf);
   }

  cat("[#.......] Reading VCF file..\n")

  readr_locale <- locale(date_names = "en", date_format = "%AD", time_format = "%AT",
                        decimal_mark = ".", grouping_mark = ",", tz = "",
                        encoding = "UTF-8", asciify = FALSE)
  
  vcf <- read_delim(filename, "\t", comment = "##", progress = FALSE, locale = readr_locale)
  vcf <- as.data.frame(vcf)
  
  vcf <- vcf[!grepl(",", vcf[ ,5]), ]
  
  if(!is.na( region_start)) vcf <- vcf[vcf[,2] >= region_start, ]
  if(!is.na( region_end)) vcf <- vcf[vcf[,2] <= region_end, ]
  
  if(is.na(maxNumberOfVariants)) maxNumberOfVariants = 1e10
  
  # Reduce number of variants
  if(!is.na(thin)) vcf <- vcf[seq(1, nrow(vcf), by = thin) ,]
  
  cat("[##......] Chromosome:  ", unique(vcf[,1]), " Mbp: ", min(vcf[,2])/1e6,
      " Region Size: ", 1e-3 * (max(vcf[,2]) - min(vcf[,2])), "kb ",
      "Num of individuals:", ncol(vcf) - 9, "\n");
  
  cat("[##......] Before filtering ", "Num of variants:", dim(vcf)[1],
      "Num of individuals:", ncol(vcf) - 9, "\n");
  
  individual_ids <- colnames(vcf)[-(1:9)]
  gt <- vcf[ , -(1:9) ]
  vcf <- vcf[ , 1:9]
  
  gt1 <- apply(gt, 2,  function(x) as.numeric(str_sub(x,1,1)))
  gt2 <- apply(gt, 2,  function(x) as.numeric(str_sub(x,3,3)))
  
  #### --- Filter by MAF ---- ####
  maf <- apply((gt1 > 0) + (gt2 > 0), 1, function(x) mean(x, na.rm = T)) / 2
  
  ok <- apply(gt1, 1, function(x) max(x, na.rm = T))
  s <- which(ok < 2)
 
  if(!is.na(min_maf)) {s <- intersect(s, which(maf >= min_maf))}
  
  if(!is.na(max_maf)) {s <- intersect(s, which(maf <= max_maf))}
  
  total_number_of_variants_within_maf <- length(s)
  
  if(length(s) > maxNumberOfVariants) s <- sort(sample(s, maxNumberOfVariants))

  gt1 <- gt1[s, ]
  gt2 <- gt2[s, ]
  gt <- gt[s, ]
  vcf <- vcf[s, ]
  maf <- maf[s]
  
  cat("[###.....] After filtering ", "Num of variants:", dim(vcf)[1],
      "Num of individuals:", ncol(gt1), "\n");

  R <- new.env()
  
  R$vcf <- vcf
  R$gt1 <- gt1
  R$gt2 <- gt2
  R$individual_ids <- individual_ids
  R$maf <- maf
  
  R$varid <- paste(R$vcf[ ,1], R$vcf[ ,2], R$vcf[ ,3], R$vcf[ ,4], R$vcf[ ,5])
  R$total_number_of_variants_within_maf <- total_number_of_variants_within_maf
  
  R
}


gene_match <- function(map_files, data_files, out_file){
  chrom <- readvcf(filename <- data_files, maxNumberOfVariants = 150, 
                       min_maf = 0, max_maf = NA)
  
  gene_map <- read.delim(map_files, header = TRUE) 
  gene_map$Chromosome <- as.numeric(sapply(gene_map$Chromosome, function(x) substring(x, 4, nchar(x))))
  
  map <- merge(gene_map, chrom$vcf, by.x = c("Chromosome", "Position.bp."),
                   by.y = c("#CHROM", "POS"))[ ,1:4]
  map$chr <- "chr"
  map <- unite(map, "Chromosome", chr, Chromosome, sep = "")
  
  write.table(map, out_file, row.names = FALSE)
}


simulation_single <- function(map_files, data_files){
  chrom <- readvcf(filename = data_files, maxNumberOfVariants = 150, 
                   min_maf = 0, max_maf = NA)
  readGeneticMap(filename = map_files)
  set.seed(1234)
  startSimulation(chrom, totalNumberOfIndividuals = 200, randomdata = 0)
  fam <- newFamilyWithOffspring("fam", 100)
  geno <- retrieveGenotypes(fam$gtindex)
  
  return(geno)
  }


simulation_mix <- function(map_files, chrom){
  readGeneticMap(filename = map_files)
  set.seed(1234)
  startSimulation(chrom, totalNumberOfIndividuals = 400, randomdata = 0)
  fam <- newFamilyWithOffspring("fam", 100)
  geno <- retrieveGenotypes(fam$gtindex)
  
  return(geno)
}


prs_single <- function (map_files, chrom_files) {
  geno <- data.frame()
  for (i in 1:length(map_files)){
    resetSimulation()
    
    chrom <- readvcf(filename = chrom_files[i], maxNumberOfVariants = 150)
    
    df <- data.frame(t(simulation_single(map_files[i], chrom_files[i])[3:102,]))
    df$chrom <- chrom$vcf[ ,1]
    df$rsids <- chrom$vcf[ ,3]
    
    geno <- rbind(geno, df)
  }
  
  return(geno)
}


prs_mix1 <- function (map_files, chrom_files1, chrom_files2) {
  geno <- data.frame()
  for (i in 1:length(map_files)){
    resetSimulation()
    
    chrom1 <- readvcf(filename = chrom_files1[i], maxNumberOfVariants = 150)
    chrom2 <- readvcf(filename = chrom_files2[i], maxNumberOfVariants = 150)
    chrom1$gt1 <- cbind(chrom1$gt1, chrom2$gt1)
    chrom1$gt2 <- cbind(chrom1$gt2, chrom2$gt2)
    chrom1$individual_ids <- c(chrom1$individual_ids, chrom2$individual_ids)
    
    df <- data.frame(t(simulation_mix(map_files[i], chrom1)[3:102,]))
    df$chrom <- chrom1$vcf[ ,1]
    df$rsids <- chrom1$vcf[ ,3]
    
    geno <- rbind(geno, df)
  }
  return(geno)
}


prs_mix2 <- function (map_files, chrom_files1, chrom_files2) {
  geno <- data.frame()
  for (i in 1:3){
    resetSimulation()
    
    chrom1 <- readvcf(filename = chrom_files1[i], maxNumberOfVariants = 150)
    chrom2 <- readvcf(filename = chrom_files2[i], maxNumberOfVariants = 150)
    chrom1$gt1 <- cbind(chrom1$gt1, chrom2$gt1)
    chrom1$gt2 <- cbind(chrom1$gt2, chrom2$gt2)
    chrom1$individual_ids <- c(chrom1$individual_ids, chrom2$individual_ids)
    
    df <- data.frame(t(simulation_mix(map_files[i], chrom1)[3:102,]))
    df$chrom <- chrom1$vcf[ ,1]
    df$rsids <- chrom1$vcf[ ,3]
    
    geno <- rbind(geno, df)
  }
  for (i in 4:(length(map_files) - 1)){
    resetSimulation()
    
    chrom1 <- readvcf(filename = chrom_files1[i + 1], maxNumberOfVariants = 150)
    chrom2 <- readvcf(filename = chrom_files2[i], maxNumberOfVariants = 150)
    chrom1$gt1 <- cbind(chrom1$gt1, chrom2$gt1)
    chrom1$gt2 <- cbind(chrom1$gt2, chrom2$gt2)
    chrom1$individual_ids <- c(chrom1$individual_ids, chrom2$individual_ids)
    
    df <- data.frame(t(simulation_mix(map_files[i + 1], chrom1)[3:102,]))
    df$chrom <- chrom1$vcf[ ,1]
    df$rsids <- chrom1$vcf[ ,3]
    
    geno <- rbind(geno, df)
  }
  return(geno)
}


prs_corr <- function (corr, geno) {
  prs_df <- merge(corr, geno, by.x = c("X.CHROM", "iID"),  by.y = c("chrom", "rsids"))
  corr_prs_df <- t(apply(prs_df[,9:length(prs_df)], 1, function(x) x[1] * x[2:length(prs_df)]))
  corr_prs <- apply(corr_prs_df, 2, sum)[1:100]
  return(corr_prs)
}


# ACB population

ACB_files <- list.files("~/Desktop/ACB_chrom", pattern = "vcf", full.names = TRUE)
ACB_map_files <- list.files("~/Desktop/ACB_map_file", pattern = "txt", full.names = TRUE)
ACB_match <- c("~/Desktop/ACB_map/ACB_map1.txt", "~/Desktop/ACB_map/ACB_map10.txt", 
               "~/Desktop/ACB_map/ACB_map11.txt", "~/Desktop/ACB_map/ACB_map12.txt", 
               "~/Desktop/ACB_map/ACB_map13.txt", "~/Desktop/ACB_map/ACB_map14.txt",
               "~/Desktop/ACB_map/ACB_map15.txt", "~/Desktop/ACB_map/ACB_map16.txt", 
               "~/Desktop/ACB_map/ACB_map17.txt", "~/Desktop/ACB_map/ACB_map18.txt", 
               "~/Desktop/ACB_map/ACB_map19.txt", "~/Desktop/ACB_map/ACB_map2.txt",
               "~/Desktop/ACB_map/ACB_map21.txt", "~/Desktop/ACB_map/ACB_map22.txt", 
               "~/Desktop/ACB_map/ACB_map3.txt", "~/Desktop/ACB_map/ACB_map4.txt", 
               "~/Desktop/ACB_map/ACB_map5.txt", "~/Desktop/ACB_map/ACB_map6.txt",  
               "~/Desktop/ACB_map/ACB_map7.txt", "~/Desktop/ACB_map/ACB_map8.txt", 
               "~/Desktop/ACB_map/ACB_map9.txt")

#for (i in 1:length(ACB_match)){
#  gene_match(ACB_map_files[i], ACB_files[i], ACB_match[i])
#  }

geno_ACB <- prs_single(ACB_match, ACB_files)


# CHS population

CHS_files <- list.files("~/Desktop/CHS_chrom", pattern = "vcf", full.names = TRUE)
CHS_map_files <- list.files("~/Desktop/CHS_map_file", pattern = "txt", full.names = TRUE)
CHS_match <- c("~/Desktop/CHS_map/CHS_map1.txt", "~/Desktop/CHS_map/CHS_map10.txt", 
               "~/Desktop/CHS_map/CHS_map11.txt", "~/Desktop/CHS_map/CHS_map13.txt", 
               "~/Desktop/CHS_map/CHS_map14.txt","~/Desktop/CHS_map/CHS_map15.txt",
               "~/Desktop/CHS_map/CHS_map16.txt", "~/Desktop/CHS_map/CHS_map17.txt",
               "~/Desktop/CHS_map/CHS_map18.txt", "~/Desktop/CHS_map/CHS_map19.txt", 
               "~/Desktop/CHS_map/CHS_map2.txt", "~/Desktop/CHS_map/CHS_map21.txt", 
               "~/Desktop/CHS_map/CHS_map22.txt", "~/Desktop/CHS_map/CHS_map3.txt", 
               "~/Desktop/CHS_map/CHS_map4.txt", "~/Desktop/CHS_map/CHS_map5.txt", 
               "~/Desktop/CHS_map/CHS_map6.txt",  "~/Desktop/CHS_map/CHS_map7.txt",
               "~/Desktop/CHS_map/CHS_map8.txt", "~/Desktop/CHS_map/CHS_map9.txt")

#for (i in 1:length(CHS_match)){
#  gene_match(CHS_map_files[i], CHS_files[i], ACB_match[i])
#  }

geno_CHS <- prs_single(CHS_match, CHS_files)


# TSI population
TSI_files <- list.files("~/Desktop/TSI_chrom", pattern = "vcf", full.names = TRUE)
TSI_map_files <- list.files("~/Desktop/TSI_map_file", pattern = "txt", full.names = TRUE)
TSI_match <- c("~/Desktop/TSI_map/TSI_map1.txt", "~/Desktop/TSI_map/TSI_map10.txt", 
               "~/Desktop/TSI_map/TSI_map11.txt", "~/Desktop/TSI_map/TSI_map12.txt", 
               "~/Desktop/TSI_map/TSI_map13.txt", "~/Desktop/TSI_map/TSI_map14.txt",
               "~/Desktop/TSI_map/TSI_map15.txt", "~/Desktop/TSI_map/TSI_map16.txt", 
               "~/Desktop/TSI_map/TSI_map17.txt", "~/Desktop/TSI_map/TSI_map18.txt", 
               "~/Desktop/TSI_map/TSI_map19.txt", "~/Desktop/TSI_map/TSI_map2.txt",
               "~/Desktop/TSI_map/TSI_map21.txt", "~/Desktop/TSI_map/TSI_map22.txt", 
               "~/Desktop/TSI_map/TSI_map3.txt", "~/Desktop/TSI_map/TSI_map4.txt", 
               "~/Desktop/TSI_map/TSI_map5.txt", "~/Desktop/TSI_map/TSI_map6.txt",  
               "~/Desktop/TSI_map/TSI_map7.txt", "~/Desktop/TSI_map/TSI_map8.txt", 
               "~/Desktop/TSI_map/TSI_map9.txt")

#for (i in 1:length(TSI_match)){
#  gene_match(TSI_map_files[i], TSI_files[i], TSI_match[i])
#}

geno_TSI <- prs_single(TSI_match, TSI_files)

# ACB and TSI population
geno_ACB_TSI <- prs_mix1(ACB_match, ACB_files, TSI_files)

# ACB and CHS population
geno_ACB_CHS <- prs_mix2(ACB_match, ACB_files, CHS_files)

# TSI and CHS population
geno_TSI_CHS <- prs_mix2(TSI_match, ACB_files, CHS_files)

#write.table(geno_TSI, "~/Desktop/prs_offspr_TSI.txt", row.names = FALSE)
#write.table(geno_CHS, "~/Desktop/prs_offspr_CHS.txt", row.names = FALSE)
#write.table(geno_ACB, "~/Desktop/prs_offspr_ACB.txt", row.names = FALSE)
#write.table(geno_TSI_CHS, "~/Desktop/prs_offspr_TSI_CHS.txt", row.names = FALSE)
#write.table(geno_ACB_CHS, "~/Desktop/prs_offspr_ACB_CHS.txt", row.names = FALSE)
#write.table(geno_TSI_ACB, "~/Desktop/prs_offspr_TSI_ACB.txt", row.names = FALSE)


# Bayes PRS with corrected beta
corr <- read.delim("~/Desktop/corrected_beta.txt")

TSI_prs <- prs_corr(corr, geno_TSI)
ACB_prs <- prs_corr(corr, geno_ACB)
CHS_prs <- prs_corr(corr, geno_CHS)
TSI_CHS_prs <- prs_corr(corr, geno_TSI_CHS)
TSI_ACB_prs <- prs_corr(corr, geno_ACB_TSI)
ACB_CHS_prs <- prs_corr(corr, geno_ACB_TSI)

