rm(list = ls())
library(ggplot2)
library(dplyr)

NumberOfTruePositives <- function( gea_results, path, suffix, LA_thresh ){
  falsePositiveDF <- data.frame(0,0,0,0,0)
  names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC","truePositive_SNPs")

  for (rep in 8:16){
    
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

#############################
### DIRECTIONAL SELECTION
#############################

bc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA.csv")
hist(bc_10000_s0.003[bc_10000_s0.003$LA!=0,]$LA,breaks = 20)
bc_10000_s0.003_LA_thresh = 6
new_bc_10000_s0.003 <- bc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map/"
bc_10000_s0.003_suffix <- "_0.003_1.12Loci.directionalSelection.csv"
bc_10000_s0.003_FP <- NumberOfTruePositives(new_bc_10000_s0.003, bc_10000_s0.003_path, bc_10000_s0.003_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.003_FP$map <- "BC Map"
bc_10000_s0.003_FP$s_max  <- "0.003"


bc_10000_s0.003_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA.sampled.csv")
hist(bc_10000_s0.003_sampled[bc_10000_s0.003_sampled$LA!=0,]$LA,breaks = 20)
bc_10000_s0.003_LA_thresh = 1
new_bc_10000_s0.003_sampled <- bc_10000_s0.003_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.003_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/"
bc_10000_s0.003_sampled_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50.csv"
bc_10000_s0.003_sampled_FP <- NumberOfTruePositives(new_bc_10000_s0.003_sampled, bc_10000_s0.003_sampled_path, bc_10000_s0.003_sampled_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.003_sampled_FP$map <- "BC Map"
bc_10000_s0.003_sampled_FP$s_max  <- "0.003"

par(mfrow  = c(3,1))
plot(bc_10000_s0.003_sampled_FP$truePositive_WZA ~ bc_10000_s0.003_sampled_FP$top)
plot(bc_10000_s0.003_sampled_FP$truePositive_TC ~ bc_10000_s0.003_sampled_FP$top)
plot(bc_10000_s0.003_sampled_FP$truePositive_SNPs ~ bc_10000_s0.003_sampled_FP$top)



bc_10000_s0.0055 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map.WZA.csv")
hist(bc_10000_s0.0055[bc_10000_s0.0055$LA!=0,]$LA,breaks = 20)
bc_10000_s0.003_LA_thresh = 20
new_bc_10000_s0.0055 <- bc_10000_s0.0055 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.0055_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map/"
bc_10000_s0.0055_suffix <- "_0.0055_1.12Loci.directionalSelection.csv"
bc_10000_s0.0055_FP <- NumberOfTruePositives(new_bc_10000_s0.0055, bc_10000_s0.0055_path, bc_10000_s0.0055_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.0055_FP$map <- "BC Map"
bc_10000_s0.0055_FP$s_max  <- "0.0055"

bc_10000_s0.0055_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map.WZA.sampled.csv")
hist(bc_10000_s0.0055_sampled[bc_10000_s0.0055_sampled$LA!=0,]$LA,breaks = 20)
bc_10000_s0.0055_LA_thresh = 1
new_bc_10000_s0.0055_sampled <- bc_10000_s0.0055_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.0055_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map_sampled/"
bc_10000_s0.0055_sampled_suffix <- "_0.0055_1.12Loci.directionalSelection_d40_n50.csv"
bc_10000_s0.0055_sampled_FP <- NumberOfTruePositives(new_bc_10000_s0.0055_sampled, bc_10000_s0.0055_sampled_path, bc_10000_s0.0055_sampled_suffix, bc_10000_s0.0055_LA_thresh)
bc_10000_s0.0055_sampled_FP$map <- "BC Map"
bc_10000_s0.0055_sampled_FP$s_max  <- "0.0055"

##################################
##################################

cline_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline.WZA.csv")
hist(cline_10000_s0.003[cline_10000_s0.003$LA!=0,]$LA,breaks = 20)
cline_10000_s0.003_LA <- 20
new_cline_10000_s0.003 <- cline_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline/"
cline_10000_s0.003_suffix <- "_0.003_3.12Loci.directionalSelection.csv"
cline_10000_s0.003_FP <- NumberOfTruePositives(new_cline_10000_s0.003, cline_10000_s0.003_path, cline_10000_s0.003_suffix, cline_10000_s0.003_LA)
cline_10000_s0.003_FP$map <- "Cline"
cline_10000_s0.003_FP$s_max <- "0.003"

cline_10000_s0.003_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline.WZA.sampled.csv")
hist(cline_10000_s0.003_sampled[cline_10000_s0.003_sampled$LA!=0,]$LA,breaks = 20)
cline_10000_s0.003_sampled_LA <- 0.6
new_cline_10000_s0.003_sampled <- cline_10000_s0.003_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.003_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled/"
cline_10000_s0.003_sampled_suffix <- "_0.003_3.12Loci.directionalSelection_d40_n50.csv"
cline_10000_s0.003_sampled_FP <- NumberOfTruePositives(new_cline_10000_s0.003_sampled, cline_10000_s0.003_sampled_path, cline_10000_s0.003_sampled_suffix, cline_10000_s0.003_sampled_LA)
cline_10000_s0.003_sampled_FP$map <- "Cline"
cline_10000_s0.003_sampled_FP$s_max <- "0.003"

cline_10000_s0.0055 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/cline.WZA.csv")
hist(cline_10000_s0.0055[cline_10000_s0.0055$LA!=0,]$LA,breaks = 20)
cline_10000_s0.0055_LA <- 50
new_cline_10000_s0.0055 <- cline_10000_s0.0055 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.0055_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/cline/"
cline_10000_s0.0055_suffix <- "_0.0055_3.12Loci.directionalSelection.csv"
cline_10000_s0.0055_FP <- NumberOfTruePositives(new_cline_10000_s0.0055, cline_10000_s0.0055_path, cline_10000_s0.0055_suffix, cline_10000_s0.0055_LA)
cline_10000_s0.0055_FP$map <- "Cline"
cline_10000_s0.0055_FP$s_max <- "0.0055"

cline_10000_s0.0055_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/cline.WZA.sampled.csv")
hist(cline_10000_s0.0055_sampled[cline_10000_s0.0055_sampled$LA!=0,]$LA,breaks = 20)
cline_10000_s0.0055_sampled_LA <- 1
new_cline_10000_s0.0055_sampled <- cline_10000_s0.0055_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.0055_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/cline_sampled/"
cline_10000_s0.0055_sampled_suffix <- "_0.0055_3.12Loci.directionalSelection_d40_n50.csv"
cline_10000_s0.0055_sampled_FP <- NumberOfTruePositives(new_cline_10000_s0.0055_sampled, cline_10000_s0.0055_sampled_path, cline_10000_s0.0055_sampled_suffix, cline_10000_s0.0055_sampled_LA)
cline_10000_s0.0055_sampled_FP$map <- "Cline"
cline_10000_s0.0055_sampled_FP$s_max <- "0.0055"

##################################
##################################

trunc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc.WZA.csv")
hist(trunc_10000_s0.003[trunc_10000_s0.003$LA!=0,]$LA,breaks = 20)
trunc_10000_s0.003_LA <- 10
new_trunc_10000_s0.003 <- trunc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc/"
trunc_10000_s0.003_suffix <- "_0.003_2.12Loci.directionalSelection.csv"
trunc_10000_s0.003_FP <- NumberOfTruePositives(new_trunc_10000_s0.003, trunc_10000_s0.003_path, trunc_10000_s0.003_suffix, trunc_10000_s0.003_LA)
trunc_10000_s0.003_FP$map <- "Truncated"
trunc_10000_s0.003_FP$s_max <- "0.003"

trunc_10000_s0.003_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc.WZA.sampled.csv")
hist(trunc_10000_s0.003_sampled[trunc_10000_s0.003_sampled$LA!=0,]$LA,breaks = 20)
trunc_10000_s0.003_sampled_LA <- 0.6
new_trunc_10000_s0.003_sampled <- trunc_10000_s0.003_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_s0.003_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled/"
trunc_10000_s0.003_sampled_suffix <- "_0.003_2.12Loci.directionalSelection.csv"
trunc_10000_s0.003_sampled_FP <- NumberOfTruePositives(new_trunc_10000_s0.003_sampled, trunc_10000_s0.003_sampled_path, trunc_10000_s0.003_sampled_suffix, trunc_10000_s0.003_sampled_LA)
trunc_10000_s0.003_sampled_FP$map <- "trunc"
trunc_10000_s0.003_sampled_FP$s_max <- "0.003"

trunc_10000_s0.0055 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/trunc.WZA.csv")
hist(trunc_10000_s0.0055[trunc_10000_s0.0055$LA!=0,]$LA,breaks = 20)
trunc_10000_s0.0055_LA <- 50
new_trunc_10000_s0.0055 <- trunc_10000_s0.0055 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_s0.0055_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/trunc/"
trunc_10000_s0.0055_suffix <- "_0.0055_3.12Loci.directionalSelection.csv"
trunc_10000_s0.0055_FP <- NumberOfTruePositives(new_trunc_10000_s0.0055, trunc_10000_s0.0055_path, trunc_10000_s0.0055_suffix, trunc_10000_s0.0055_LA)
trunc_10000_s0.0055_FP$map <- "trunc"
trunc_10000_s0.0055_FP$s_max <- "0.0055"

trunc_10000_s0.0055_sampled <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/trunc.WZA.sampled.csv")
hist(trunc_10000_s0.0055_sampled[trunc_10000_s0.0055_sampled$LA!=0,]$LA,breaks = 20)
trunc_10000_s0.0055_sampled_LA <- 1
new_trunc_10000_s0.0055_sampled <- trunc_10000_s0.0055_sampled %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_s0.0055_sampled_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/trunc_sampled/"
trunc_10000_s0.0055_sampled_suffix <- "_0.0055_3.12Loci.directionalSelection_d40_n50.csv"
trunc_10000_s0.0055_sampled_FP <- NumberOfTruePositives(new_trunc_10000_s0.0055_sampled, trunc_10000_s0.0055_sampled_path, trunc_10000_s0.0055_sampled_suffix, trunc_10000_s0.0055_sampled_LA)
trunc_10000_s0.0055_sampled_FP$map <- "trunc"
trunc_10000_s0.0055_sampled_FP$s_max <- "0.0055"


sampled <- rbind(bc_10000_s0.003_sampled_FP,bc_10000_s0.0055_sampled_FP,cline_10000_s0.003_sampled_FP,cline_10000_s0.0055_sampled_FP)
sampled$s_max <- factor(sampled$s_max, levels = c("0.003", "0.0055"), labels = c(expression(italic(s[Max])*" = 0.003"), expression(italic(s[Max])*" = 0.0055")))
sampled$map <- factor(sampled$map, levels = c("BC Map", "Cline"), labels = c(expression(italic("BC Map")), expression(italic("Cline"))))

all_data <- rbind(bc_10000_s0.003_FP,bc_10000_s0.0055_FP,cline_10000_s0.003_FP,cline_10000_s0.0055_FP)
all_data$s_max <- factor(all_data$s_max, levels = c("0.003", "0.0055"), labels = c(expression(italic(s[Max])*" = 0.003"), expression(italic(s[Max])*" = 0.0055")))
all_data$map <- factor(all_data$map, levels = c("BC Map", "Cline"), labels = c(expression(italic("BC Map")), expression(italic("Cline"))))

sampled_data_plot <- ggplot(data = sampled, aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = sampled, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = sampled, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = sampled, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  facet_grid(map~s_max, labeller = label_parsed)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )



all_data_plot <- ggplot(data = all_data, aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = all_data, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = all_data, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = all_data, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  facet_grid(map~s_max, labeller = label_parsed)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )

png("~/work/GEA/simulations/directionalSelection/PowerPlot_allData.png", width = 8, height = 6, units = "in", res = 300)
print(all_data_plot)
dev.off()

png("~/work/GEA/simulations/directionalSelection/PowerPlot_sampledData.png", width = 8, height = 6, units = "in", res = 300)
print(sampled_data_plot)
dev.off()






all_data_plot <- ggplot(data = trunc_10000_s0.003_sampled_FP, aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = trunc_10000_s0.003_sampled_FP, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = trunc_10000_s0.003_sampled_FP, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = trunc_10000_s0.003_sampled_FP, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  facet_grid(map~s_max, labeller = label_parsed)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )











bc_10000_s0.003_sampledRandom <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.sampled_randomPops.csv")
hist(bc_10000_s0.003_sampledRandom[bc_10000_s0.003_sampledRandomRandom$LA!=0,]$LA,breaks = 20)
bc_10000_s0.003_LA_thresh = 1
new_bc_10000_s0.003_sampledRandom <- bc_10000_s0.003_sampledRandom %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.003_sampledRandom_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/"
bc_10000_s0.003_sampledRandom_suffix <- "_0.003_1.12Loci.directionalSelection_d40_n50.csv"
bc_10000_s0.003_sampledRandom_FP <- NumberOfTruePositives(new_bc_10000_s0.003_sampledRandom, bc_10000_s0.003_sampledRandom_path, bc_10000_s0.003_sampledRandom_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.003_sampledRandom_FP$map <- "BC Map Random Pops"
bc_10000_s0.003_sampledRandom_FP$s_max  <- "0.003"


new <- rbind(bc_10000_s0.003_sampled_FP, bc_10000_s0.003_sampledRandom_FP)

sampled_data_plot <- ggplot(data = rbind(bc_10000_s0.003_sampled_FP, bc_10000_s0.003_sampledRandom_FP), aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = new, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = new, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = new, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  facet_grid(map~s_max)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )





##################################
##    STABILISNG SELECTION
##################################



ggplot(data = bc_10000_Vs192_sampledRandom_FP, aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = bc_10000_Vs192_sampledRandom_FP, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = bc_10000_Vs192_sampledRandom_FP, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = bc_10000_Vs192_sampledRandom_FP, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )

##################################333

ggplot(data = bc_10000_Vs192_sampledRandom, aes( x = position, y = X_Z, col = LA > 0.01))+
  geom_point()+
  facet_wrap(~rep)

ggplot(data = bc_10000_Vs192_sampledRandom, aes( x = position, y = -log10(top_candidate_p), col = LA > 0.01))+
  geom_point()+
  facet_wrap(~rep)

ggplot(data = bc_10000_Vs192_sampledRandom, aes( x  = X_Z, fill = LA > 0.01))+
  geom_density()  
