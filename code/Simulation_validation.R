# simulation for predict performance

snp <- c(0:2)
set.seed(123)
#set 3 types of SNPs: maf = 0.1, 0.2, 0.4
p=c(0.1,0.2,0.4)
q=1-p
genotype_freq=cbind(q^2,2*p*q,p^2)
#n=500
genotypes1 <- matrix(c(sample(snp,2000000,replace = T, prob=genotype_freq[1,]),
                      sample(snp,2000000,replace = T, prob=genotype_freq[2,]),
                      sample(snp,1000000,replace = T, prob=genotype_freq[3,])),ncol = 10000, nrow = 500)

pop1 <- data.frame(cbind(c(1:500),sample(c(0,1),500,replace = T, prob=c(0.5,0.5)),genotypes1))

set.seed(1234)
#n=1000
genotypes2 <- matrix(c(sample(snp,4000000,replace = T, prob=genotype_freq[1,]),
                      sample(snp,4000000,replace = T, prob=genotype_freq[2,]),
                      sample(snp,2000000,replace = T, prob=genotype_freq[3,])),ncol = 10000, nrow = 1000)

pop1 <- data.frame(cbind(c(1:1000),sample(c(0,1),1000,replace = T, prob=c(0.5,0.5)),genotypes2))


