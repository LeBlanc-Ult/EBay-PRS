setwd("/Users/Lydia/UCSD/SP21/CSE284/Proj")
ACB <- read.table(file = "prs_offspr_ACB.txt", sep = " ", header = T)
CHS <- read.table(file = "prs_offspr_CHS.txt", sep = " ", header = T)
TSI <- read.table(file = "prs_offspr_TSI.txt", sep = " ", header = T)
TSI_ACB <- read.table(file = "prs_offspr_TSI_ACB.txt", sep = " ", header = T)
TSI_CHS <- read.table(file = "prs_offspr_TSI_CHS.txt", sep = " ", header = T)
ACB_CHS <- read.table(file = "prs_offspr_ACB_CHS.txt", sep = " ", header = T)

library(ggplot2)
library(ggpubr)

p1 <- ggplot(ACB, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),      
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") + 
  xlab("PRS of ACB x ACB") +
  xlim(-0.1,0.1)


p2 <- ggplot(CHS, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),      
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") + 
  xlab("PRS of CHS x CHS")+
  xlim(-0.1,0.1)

p3 <- ggplot(TSI, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),      
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +  
  xlab("PRS of TSI x TSI")+
  xlim(-0.1,0.1)



p4 <- ggplot(TSI_ACB, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),      
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") + 
  xlab("PRS of TSI x ACB") +
  xlim(-0.1,0.1)

p5 <- ggplot(TSI_CHS, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),      
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +  
  xlab("PRS of TSI x CHS") +
  xlim(-0.1,0.1)


p6 <- ggplot(ACB_CHS, aes(x=prs)) + 
  geom_histogram(aes(y=..density..),    
                 binwidth=.005,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") + 
  xlab("PRS of ACB x CHS") +
  xlim(-0.1,0.1)
ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2,labels=c("A","B","C","D","E","F"))

quantile(ACB$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975))
round(quantile(CHS$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975)),5)
quantile(TSI$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975))
round(quantile(TSI_ACB$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975)),5)
round(quantile(TSI_CHS$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975)),5)
round(quantile(ACB_CHS$prs, probs = c(0.025, 0.2, 0.5, 0.8, 0.975)),5)
