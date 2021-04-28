rm(list=ls())


#gea_results <- new_bc_10_s0.003
#LA_thresh =1 

#path <- "~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10/"
#suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50_i10.csv"

NumberOfTruePositives <- function( gea_results, path, suffix, LA_thresh ){
  falsePositiveDF <- data.frame(0,0,0,0,0)
  names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC","truePositive_SNPs")
  
  for (rep in 1:20){
    
    temp <- gea_results[gea_results$rep == rep,]
    numberLA <- length(droplevels(temp[temp$LA > LA_thresh,]$gene))
    LA_genes <- droplevels(unique(temp[temp$LA > LA_thresh,]$gene))
    
    snps <- read.csv(paste(path, rep, suffix, sep = ''))
    
    snps <- snps[ with(snps, order(ave(pop_k_tau_p_value, gene, FUN = min), pop_k_tau_p_value)),]
    snp_based_genes <- unique( snps$gene )
    
    for (n in 1:50){
      
      WZA_slice <- temp[temp$Z_empirical_p <= (n+1)/1000,]
      TC_slice <- temp[temp$TC_empirical_p <= (n+1)/1000,]
      SNP_slice <- snp_based_genes[1:n]
      
      #      SNP_slice <- snps[snps$rank <= n,] ## Grab the top n SNPs from the genome scan
      
      #      numberLA_SNP_slice <- length(unname(unique(SNP_slice[-log10(SNP_slice$LA) > 3.0,]$gene)))
      
      truePositive_WZA = sum( WZA_slice$LA > LA_thresh )/numberLA
      truePositive_TC = sum( TC_slice$LA > LA_thresh )/numberLA
      truePositive_SNPs = sum(SNP_slice%in%LA_genes)/numberLA
      
      falsePositiveDF[(50*(rep-1))+n,] <- c(rep, n, truePositive_WZA, truePositive_TC, truePositive_SNPs)
      
    }
  }
  falsePositiveDF
}

b10 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10.WZA.csv"))
b10$bLen <-"10"
#b10$top_candidate_p<- correctTC(b10)

b100 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_100.WZA.csv"))
b100$bLen <-"100"
#b100$top_candidate_p<- correctTC(b100)

b1000 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_1000.WZA.csv"))
b1000$bLen <-"1000"
#b1000$top_candidate_p<- correctTC(b1000)


b10000 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA.sampled.csv"))
b10000$bLen <-"10000"
#b10000$top_candidate_p<- correctTC(b10000)

b10 <- b10[order(b10$rep,b10$gene),]
b100 <- b100[order(b100$rep,b100$gene),]
b1000 <- b1000[order(b1000$rep,b1000$gene),]
b10000 <- b10000[order(b10000$rep,b10000$gene),]

blocks <- rbind(b10, b100, b1000, b10000)
blocks <- blocks[blocks$rep <=2,]


## Make a density plot showing recombination distances and the variance in TC and Z
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(reshape2)



pal <- wes_palette("FantasticFox1")[c(1,3,4,5)]
pal <- c("#d7191c","#fdae61","#abd9e9","#2c7bb6")
blocks_trunc <- blocks[blocks$position > 8e6,]

blocksDensity_Z <- ggplot(data = blocks_trunc, aes(x = Z, fill = bLen))+
  geom_density(alpha = 0.6)+
  theme_bw()+
  scale_fill_manual("Block Length",values= pal)+
  scale_x_continuous(expression(Z[W]))


blocksDensity_TC <- ggplot(data = blocks_trunc, aes(x = -log10(top_candidate_p), fill = bLen))+
  geom_density(alpha = 0.6)+
  theme_bw()+
  scale_fill_manual("Block Length",values= pal)+
  scale_x_continuous("Top Candidate Index",limits = c(0,5))

png("~/work/GEA/simulations/directionalSelection/StatRecombinationComparisons_sampled.png", width = 10, height = 4, units = "in", res = 300)
ggarrange(blocksDensity_Z, blocksDensity_TC,common.legend = T, legend = "right")
dev.off()



WZA <- ggplot(data = blocks, aes(x = position/1e6, y = Z,col =LA>1))+
  geom_vline(xintercept = c(2,4,6,8), lty = 2, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("black","red"))+
  scale_x_continuous("Position")+
  facet_grid(bLen~rep, scales = "free_y")+
  theme_bw()
TC <- ggplot(data = blocks, aes(x = position/1e6, y = -log10(top_candidate_p),col =LA>1))+
  geom_vline(xintercept = c(2,4,6,8), lty = 2, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("black","red"))+
  scale_x_continuous("Position")+
  facet_grid(bLen~rep, scales = "free_y")+
  theme_bw()



r_1_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10/1_0.003_1.12Loci.directionalSelection_d40_n50_i10.csv")  
r_1_10$bLen <-"10"
r_1_10$rep <- "1"

r_1_100 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_100/1_0.003_1.12Loci.directionalSelection_d40_n50_i100.csv")  
r_1_100$bLen <-"100"
r_1_100$rep <- "1"

r_1_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_1000/1_0.003_1.12Loci.directionalSelection_d40_n50_i1000.csv")  
r_1_1000$bLen <-"1000"
r_1_1000$rep <- "1"

r_1_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/1_0.003_1.12Loci.directionalSelection_d40_n50.csv")  
r_1_10000$bLen <-"10000"
r_1_10000$rep <- "1"

r_2_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10/2_0.003_1.12Loci.directionalSelection_d40_n50_i10.csv")  
r_2_10$bLen <-"10"
r_2_10$rep <- "2"

r_2_100 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_100/2_0.003_1.12Loci.directionalSelection_d40_n50_i100.csv")  
r_2_100$bLen <-"100"
r_2_100$rep <- "2"

r_2_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_1000/2_0.003_1.12Loci.directionalSelection_d40_n50_i1000.csv")  
r_2_1000$bLen <-"1000"
r_2_1000$rep <- "2"

r_2_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/2_0.003_1.12Loci.directionalSelection_d40_n50.csv")  
r_2_10000$bLen <-"10000"
r_2_10000$rep <- "2"

SNPs <- ggplot(data = rbind(r_1_10, r_1_100,r_1_1000, r_1_10000, r_2_10,r_2_100, r_2_1000, r_2_10000), aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA > 1))+
  geom_vline(xintercept = c(2,4,6,8), lty = 2, alpha = 0.3)+
  geom_point(alpha = 0.5)+
  theme_bw()+
  scale_x_continuous("Position")+
  scale_color_manual(values = c("black","red"))+
  facet_grid(bLen~rep, scales = "free_y")

png("~/work/GEA/simulations/directionalSelection/StatRecombinationComparisons_AcrossReps_sampled.png", width = 20, height = 8, units = "in", res = 300)
ggarrange(SNPs, TC, WZA, legend = "right", nrow =1, ncol = 3, labels = "AUTO")
dev.off()




######################
######################



bc_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA.sampled.csv")
#bc_10000$top_candidate_p<- correctTC(bc_10000)
hist(bc_10000[bc_10000$LA!=0,]$LA,breaks = 20)
bc_10000_s0.003_LA_thresh = 1
new_bc_10000_s0.003 <- bc_10000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/"
bc_10000_s0.003_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50.csv"
bc_10000_s0.003_FP <- NumberOfTruePositives(new_bc_10000_s0.003, bc_10000_s0.003_path, bc_10000_s0.003_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.003_FP$map <- "BC Map"
bc_10000_s0.003_FP$s_max  <- "0.003"
bc_10000_s0.003_FP$blockSize <- "10000"
bc_10000_s0.003_FP_melt <- melt(bc_10000_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))





bc_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_1000.WZA.csv")
#bc_1000$top_candidate_p<- correctTC(bc_1000)
hist(bc_1000[bc_1000$LA!=0,]$LA,breaks = 20)
bc_1000_s0.003_LA_thresh = 1
new_bc_1000_s0.003 <- bc_1000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_1000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_1000/"
bc_1000_s0.003_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50_i1000.csv"
bc_1000_s0.003_FP <- NumberOfTruePositives(new_bc_1000_s0.003, bc_1000_s0.003_path, bc_1000_s0.003_suffix, bc_1000_s0.003_LA_thresh)
bc_1000_s0.003_FP$map <- "BC Map"
bc_1000_s0.003_FP$s_max  <- "0.003"
bc_1000_s0.003_FP$blockSize <- "1000"
bc_1000_s0.003_FP_melt <- melt(bc_1000_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))


bc_100 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_100.WZA.csv")
#bc_100$top_candidate_p<- correctTC(bc_100)
hist(bc_100[bc_100$LA!=0,]$LA,breaks = 20)
bc_100_s0.003_LA_thresh = 1
new_bc_100_s0.003 <- bc_100 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_100_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_100/"
bc_100_s0.003_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50_i100.csv"
bc_100_s0.003_FP <- NumberOfTruePositives(new_bc_100_s0.003, bc_100_s0.003_path, bc_100_s0.003_suffix, bc_100_s0.003_LA_thresh)
bc_100_s0.003_FP$map <- "BC Map"
bc_100_s0.003_FP$s_max  <- "0.003"
bc_100_s0.003_FP$blockSize <- "100"
bc_100_s0.003_FP_melt <- melt(bc_100_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))




bc_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10.WZA.csv")
#bc_10$top_candidate_p<- correctTC(bc_10)
hist(bc_10[bc_10$LA!=0,]$LA,breaks = 20)
bc_10_s0.003_LA_thresh = 1
new_bc_10_s0.003 <- bc_10 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_sampled_10/"
bc_10_s0.003_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50_i10.csv"
bc_10_s0.003_FP <- NumberOfTruePositives(new_bc_10_s0.003, bc_10_s0.003_path, bc_10_s0.003_suffix, bc_10_s0.003_LA_thresh)
bc_10_s0.003_FP$map <- "BC Map"
bc_10_s0.003_FP$s_max  <- "0.003"
bc_10_s0.003_FP$blockSize <- "10"
names(bc_10_s0.003_FP)
bc_10_s0.003_FP_melt <- melt(bc_10_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs" ))

bc_blockers <- rbind(bc_10_s0.003_FP_melt,bc_100_s0.003_FP_melt,bc_1000_s0.003_FP_melt,bc_10000_s0.003_FP_melt)



bcResult <- ggplot(data = bc_blockers, aes(x = top, y = value))+
  stat_summary(data = bc_blockers, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = bc_blockers, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = bc_blockers, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
#  stat_summary(data = blockers, aes(y = value, x = top, col = blockSize), fun.y = mean, geom = "line", lwd = 2)+
  facet_grid(value~.)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  scale_x_continuous("Top # Genes")+
  theme_bw()


# 





######################
######################



cline_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline.WZA.sampled.csv")
hist(cline_10000[cline_10000$LA!=0,]$LA,breaks = 20)
cline_10000_s0.003_LA_thresh = 1
new_cline_10000_s0.003 <- cline_10000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled/"
cline_10000_s0.003_suffix <- "_0.003_3.12Loci.directionalSelection_d40_n50.csv"
cline_10000_s0.003_FP <- NumberOfTruePositives(new_cline_10000_s0.003, cline_10000_s0.003_path, cline_10000_s0.003_suffix, cline_10000_s0.003_LA_thresh)
cline_10000_s0.003_FP$map <- "cline Map"
cline_10000_s0.003_FP$s_max  <- "0.003"
cline_10000_s0.003_FP$blockSize <- "10000"
cline_10000_s0.003_FP_melt <- melt(cline_10000_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))





cline_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_1000.WZA.csv")
hist(cline_1000[cline_1000$LA!=0,]$LA,breaks = 20)
cline_1000_s0.003_LA_thresh = 1
new_cline_1000_s0.003 <- cline_1000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_1000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_1000/"
cline_1000_s0.003_suffix <- "_0.003_3.12Loci.directionalSelection_d40_n50_i1000.csv"
cline_1000_s0.003_FP <- NumberOfTruePositives(new_cline_1000_s0.003, cline_1000_s0.003_path, cline_1000_s0.003_suffix, cline_1000_s0.003_LA_thresh)
cline_1000_s0.003_FP$map <- "cline Map"
cline_1000_s0.003_FP$s_max  <- "0.003"
cline_1000_s0.003_FP$blockSize <- "1000"
cline_1000_s0.003_FP_melt <- melt(cline_1000_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))


cline_100 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_100.WZA.csv")
hist(cline_100[cline_100$LA!=0,]$LA,breaks = 20)
cline_100_s0.003_LA_thresh = 1
new_cline_100_s0.003 <- cline_100 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_100_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_100/"
cline_100_s0.003_suffix <- "_0.003_3.12Loci.directionalSelection_d40_n50_i100.csv"
cline_100_s0.003_FP <- NumberOfTruePositives(new_cline_100_s0.003, cline_100_s0.003_path, cline_100_s0.003_suffix, cline_100_s0.003_LA_thresh)
cline_100_s0.003_FP$map <- "cline Map"
cline_100_s0.003_FP$s_max  <- "0.003"
cline_100_s0.003_FP$blockSize <- "100"
cline_100_s0.003_FP_melt <- melt(cline_100_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs"))




cline_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_10.WZA.csv")
hist(cline_10[cline_10$LA!=0,]$LA,breaks = 20)
cline_10_s0.003_LA_thresh = 1
new_cline_10_s0.003 <- cline_10 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/cline_recombinationVariation/cline_0.003_sampled_10/"
cline_10_s0.003_suffix <- "_0.003_3.12Loci.directionalSelection_d40_n50_i10.csv"
cline_10_s0.003_FP <- NumberOfTruePositives(new_cline_10_s0.003, cline_10_s0.003_path, cline_10_s0.003_suffix, cline_10_s0.003_LA_thresh)
cline_10_s0.003_FP$map <- "cline Map"
cline_10_s0.003_FP$s_max  <- "0.003"
cline_10_s0.003_FP$blockSize <- "10"
names(cline_10_s0.003_FP)
cline_10_s0.003_FP_melt <- melt(cline_10_s0.003_FP, id = c("rep","top","map","s_max", "truePositive_WZA","truePositive_TC","truePositive_SNPs" ))

cline_blockers <- rbind(cline_10_s0.003_FP_melt,cline_100_s0.003_FP_melt,cline_1000_s0.003_FP_melt,cline_10000_s0.003_FP_melt)



clineResult <- ggplot(data = cline_blockers, aes(x = top, y = value))+
  stat_summary(data =cline_blockers, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = cline_blockers, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = cline_blockers, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  #  stat_summary(data = blockers, aes(y = value, x = top, col = blockSize), fun.y = mean, geom = "line", lwd = 2)+
  facet_grid(value~.)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  scale_x_continuous("Top # Genes")+
  theme_bw()


blockPlot <- ggarrange(blocksDensity_Z, blocksDensity_TC,common.legend = T, legend = "right", nrow= 2, ncol = 1,  labels = c("A","B"))

linePlot <-  ggarrange( bcResult, clineResult, nrow= 1, ncol=2,common.legend = T, labels = c("C","D"))


ggarrange(blockPlot,linePlot)


# r_1_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_10_sampled/1_0.003_1.12Loci.directionalSelection_10bp_p40_n50.csv")  
# r_1_10$bLen <-"10"
# r_1_10$rep <- "1"
# 
# r_1_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_1000_sampled/1_0.003_1.12Loci.directionalSelection_1000bp_p40_n50.csv")  
# r_1_1000$bLen <-"1000"
# r_1_1000$rep <- "1"
# 
# r_1_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_10000_sampled/1_0.003_1.12Loci.directionalSelection_10000bp_p40_n50.csv")  
# r_1_10000$bLen <-"10000"
# r_1_10000$rep <- "1"
# 
# r_2_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_10_sampled/2_0.003_1.12Loci.directionalSelection_10bp_p40_n50.csv")  
# r_2_10$bLen <-"10"
# r_2_10$rep <- "2"
# 
# r_2_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_1000_sampled/2_0.003_1.12Loci.directionalSelection_1000bp_p40_n50.csv")  
# r_2_1000$bLen <-"1000"
# r_2_1000$rep <- "2"
# 
# r_2_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/test/BC_Map_0.003_10000_sampled/2_0.003_1.12Loci.directionalSelection_10000bp_p40_n50.csv")  
# r_2_10000$bLen <-"10000"
# r_2_10000$rep <- "2"
# 
# SNPs <- ggplot(data = rbind(r_1_10, r_1_1000, r_1_10000, r_2_10, r_2_1000, r_2_10000), aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA < 1))+
#   geom_vline(xintercept = c(2,4,6,8), lty = 2, alpha = 0.3)+
#   geom_point(alpha = 0.5)+
#   theme_bw()+
#   scale_x_continuous("Position")+
#   scale_color_manual(values = c("black","red"))+
#   facet_grid(bLen~rep, scales = "free_y")
# 
# library(ggpubr)
# 
# png("~/work/GEA/simulations/directionalSelection/BC_map_comparison_sampled.png", width = 12, height = 6, units = "in", res = 300)
# ggarrange( SNPs, WZA, TC ,nrow= 1, ncol=3,common.legend = T)
# dev.off()
# 
# 


