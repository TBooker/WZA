rm(list = ls ())

library(ggplot2)
library(dplyr)

NumberOfTruePositives <- function( gea_results, path, suffix, LA_thresh ){
  falsePositiveDF <- data.frame(0,0,0,0,0)
  names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC","truePositive_SNPs")
  
  for (rep in 1:20){
    
    temp <- gea_results[gea_results$rep == rep,]
    temp <- temp[temp$LA != -99,]
    numberLA <- length(droplevels(temp[temp$LA > LA_thresh,]$gene))
    LA_genes <- droplevels(unique(temp[temp$LA > LA_thresh,]$gene))
    
    snps <- read.csv(paste(path, rep, suffix, sep = ''))
    snps <- snps[snps$LA != -99,]
    
    snps <- snps[ with(snps, order(ave(geno_k_tau_p_value, gene, FUN = min), geno_k_tau_p_value)),]
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


summarise <- function(top_hits){  
  wza <-  as.data.frame(aggregate( truePositive_WZA ~ top, data = top_hits, 
                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)))) ) 
  wza[,2:3] = unlist(wza$truePositive_WZA)
  names( wza ) <- c("top", "mean", "se")
  wza$test <- "WZA"
  tc <-  aggregate( truePositive_TC ~ top, data = top_hits,
                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
  tc[,2:3] = unlist(tc$truePositive_TC)
  names( tc ) <- c("top", "mean", "se")
  tc$test <- "Top-Candidate"
  snp <- aggregate( truePositive_SNPs ~ top, data = top_hits,
                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
  snp[,2:3] = unlist(snp$truePositive_SNPs)
  names( snp ) <- c("top", "mean", "se")
  snp$test <- "SNP-based"
  
  return( rbind( wza, tc, snp) )
}

############################
##    BC Map

bc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA.csv")
hist(bc_10000_s0.003[bc_10000_s0.003$LA>0,]$LA, breaks = 20)
bc_10000_s0.003_LA_thresh = 0.005
new_bc_10000_s0.003 <- bc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map/"
bc_10000_s0.003_suffix <- "_0.003_directional_i10000.csv"
bc_10000_s0.003_FP <- NumberOfTruePositives(new_bc_10000_s0.003, bc_10000_s0.003_path, bc_10000_s0.003_suffix, bc_10000_s0.003_LA_thresh)
bc_10000_s0.003_summary <- summarise( bc_10000_s0.003_FP )
bc_10000_s0.003_summary$map <- "BC Map"
bc_10000_s0.003_summary$selection <- "Directional"

bc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/BC_Map.WZA.csv")
hist(bc_10000_Vs192[bc_10000_Vs192$LA>0,]$LA,breaks = 20)
bc_10000_Vs192_LA_thresh = 0.005
new_bc_10000_Vs192 <- bc_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
bc_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/BC_Map/"
bc_10000_Vs192_suffix <- "_0.5_192_d40n50_i10000.csv"
bc_10000_Vs192_FP <- NumberOfTruePositives(new_bc_10000_Vs192, bc_10000_Vs192_path, bc_10000_Vs192_suffix, bc_10000_Vs192_LA_thresh)
bc_10000_Vs192_FP$map <- "BC Map"
bc_10000_Vs192_FP$s_max  <- "0.003"
bc_10000_Vs192_summary <- summarise( bc_10000_Vs192_FP )
bc_10000_Vs192_summary$map <- "BC Map"
bc_10000_Vs192_summary$selection <- "Stabilising"

########################## Cline


cline_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline.WZA.csv")
hist(cline_10000_s0.003[cline_10000_s0.003$LA>0,]$LA, breaks = 20)
cline_10000_s0.003_LA_thresh = 0.005
new_cline_10000_s0.003 <- cline_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline/"
cline_10000_s0.003_suffix <- "_0.003_directional_i10000.csv"
cline_10000_s0.003_FP <- NumberOfTruePositives(new_cline_10000_s0.003, cline_10000_s0.003_path, cline_10000_s0.003_suffix, cline_10000_s0.003_LA_thresh)
cline_10000_s0.003_FP$map <- "cline"
cline_10000_s0.003_FP$s_max  <- "0.003"
cline_10000_s0.003_summary <- summarise( cline_10000_s0.003_FP )
cline_10000_s0.003_summary$map <- "Cline"
cline_10000_s0.003_summary$selection <- "Directional"


cline_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/cline.WZA.csv")
hist(cline_10000_Vs192[cline_10000_Vs192$LA>0,]$LA,breaks = 20)
cline_10000_Vs192_LA_thresh = 0.005
new_cline_10000_Vs192 <- cline_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
cline_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/cline/"
cline_10000_Vs192_suffix <- "_0.5_192_i10000.csv"
cline_10000_Vs192_FP <- NumberOfTruePositives(new_cline_10000_Vs192, cline_10000_Vs192_path, cline_10000_Vs192_suffix, cline_10000_Vs192_LA_thresh)
cline_10000_Vs192_FP$map <- "cline"
cline_10000_Vs192_FP$s_max  <- "0.003"
cline_10000_Vs192_summary <- summarise( cline_10000_Vs192_FP )
cline_10000_Vs192_summary$map <- "Cline"
cline_10000_Vs192_summary$selection <- "Stabilising"



################## Trunc


trunc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc.WZA.csv")
hist(trunc_10000_s0.003[trunc_10000_s0.003$LA>0,]$LA, breaks = 20)
trunc_10000_s0.003_LA_thresh = 0.005
new_trunc_10000_s0.003 <- trunc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc/"
trunc_10000_s0.003_suffix <- "_0.003_directional_i10000.csv"
trunc_10000_s0.003_FP <- NumberOfTruePositives(new_trunc_10000_s0.003, trunc_10000_s0.003_path, trunc_10000_s0.003_suffix, trunc_10000_s0.003_LA_thresh)
trunc_10000_s0.003_FP$map <- "trunc"
trunc_10000_s0.003_FP$s_max  <- "0.003"
trunc_10000_s0.003_summary <- summarise( trunc_10000_s0.003_FP )
trunc_10000_s0.003_summary$map <- "Truncated"
trunc_10000_s0.003_summary$selection <- "Directional"



trunc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/trunc.WZA.csv")
hist(trunc_10000_Vs192[trunc_10000_Vs192$LA>0,]$LA,breaks = 20)
trunc_10000_Vs192_LA_thresh = 0.005
new_trunc_10000_Vs192 <- trunc_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
trunc_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/trunc/"
trunc_10000_Vs192_suffix <- "_0.5_192_i10000.csv"
trunc_10000_Vs192_FP <- NumberOfTruePositives(new_trunc_10000_Vs192, trunc_10000_Vs192_path, trunc_10000_Vs192_suffix, trunc_10000_Vs192_LA_thresh)
trunc_10000_Vs192_FP$map <- "trunc"
trunc_10000_Vs192_FP$s_max  <- "0.003"
trunc_10000_Vs192_summary <- summarise( trunc_10000_Vs192_FP )
trunc_10000_Vs192_summary$map <- "Truncated"
trunc_10000_Vs192_summary$selection <- "Stabilising"




###############################
##  Sampled BC


sampled_bc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled.WZA.csv")
hist(sampled_bc_10000_s0.003[sampled_bc_10000_s0.003$LA>0,]$LA, breaks = 20)
sampled_bc_10000_s0.003_LA_thresh = 0.005
new_sampled_bc_10000_s0.003 <- sampled_bc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_bc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/"
sampled_bc_10000_s0.003_suffix <- "_0.003_directional_d40n50_i10000.csv"
sampled_bc_10000_s0.003_FP <- NumberOfTruePositives(new_sampled_bc_10000_s0.003, sampled_bc_10000_s0.003_path, sampled_bc_10000_s0.003_suffix, sampled_bc_10000_s0.003_LA_thresh)
sampled_bc_10000_s0.003_FP$map <- "BC Map"
sampled_bc_10000_s0.003_FP$s_max  <- "0.003"
sampled_bc_10000_s0.003_summary <- summarise( sampled_bc_10000_s0.003_FP )
sampled_bc_10000_s0.003_summary$map <- "BC Map"
sampled_bc_10000_s0.003_summary$selection <- "Directional"

sampled_bc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/BC_Map_sampled.WZA.csv")
hist(sampled_bc_10000_Vs192[sampled_bc_10000_Vs192$LA>0,]$LA, breaks = 20)
sampled_bc_10000_Vs192_LA_thresh = 0.005
new_sampled_bc_10000_Vs192 <- sampled_bc_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_bc_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/BC_Map_sampled/"
sampled_bc_10000_Vs192_suffix <- "_0.5_192_d40n50_i10000.csv"
sampled_bc_10000_Vs192_FP <- NumberOfTruePositives(new_sampled_bc_10000_Vs192, sampled_bc_10000_Vs192_path, sampled_bc_10000_Vs192_suffix, sampled_bc_10000_Vs192_LA_thresh)
sampled_bc_10000_Vs192_FP$map <- "BC Map"
sampled_bc_10000_Vs192_FP$s_max  <- "0.003"
sampled_bc_10000_Vs192_summary <- summarise( sampled_bc_10000_Vs192_FP )
sampled_bc_10000_Vs192_summary$map <- "BC Map"
sampled_bc_10000_Vs192_summary$selection <- "Stabilising"


###############################
##  Sampled cline


sampled_cline_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled.WZA.csv")
hist(sampled_cline_10000_s0.003[sampled_cline_10000_s0.003$LA>0,]$LA, breaks = 20)
sampled_cline_10000_s0.003_LA_thresh = 0.005
new_sampled_cline_10000_s0.003 <- sampled_cline_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_cline_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled/"
sampled_cline_10000_s0.003_suffix <- "_0.003_directional_d40n50_i10000.csv"
sampled_cline_10000_s0.003_FP <- NumberOfTruePositives(new_sampled_cline_10000_s0.003, sampled_cline_10000_s0.003_path, sampled_cline_10000_s0.003_suffix, sampled_cline_10000_s0.003_LA_thresh)
sampled_cline_10000_s0.003_FP$map <- "cline Map"
sampled_cline_10000_s0.003_FP$s_max  <- "0.003"
sampled_cline_10000_s0.003_summary <- summarise( sampled_cline_10000_s0.003_FP )
sampled_cline_10000_s0.003_summary$map <- "Cline"
sampled_cline_10000_s0.003_summary$selection <- "Directional"

sampled_cline_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/cline_sampled_1000.WZA.csv")
hist(sampled_cline_10000_Vs192[sampled_cline_10000_Vs192$LA>0,]$LA, breaks = 20)
sampled_cline_10000_Vs192_LA_thresh = 0.005
new_sampled_cline_10000_Vs192 <- sampled_cline_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_cline_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/cline_sampled_1000/"
sampled_cline_10000_Vs192_suffix <- "_0.5_192_d40n50_i1000.csv"
sampled_cline_10000_Vs192_FP <- NumberOfTruePositives(new_sampled_cline_10000_Vs192, sampled_cline_10000_Vs192_path, sampled_cline_10000_Vs192_suffix, sampled_cline_10000_Vs192_LA_thresh)
sampled_cline_10000_Vs192_FP$map <- "cline Map"
sampled_cline_10000_Vs192_FP$s_max  <- "0.003"
sampled_cline_10000_Vs192_summary <- summarise( sampled_cline_10000_Vs192_FP )
sampled_cline_10000_Vs192_summary$map <- "Cline"
sampled_cline_10000_Vs192_summary$selection <- "Stabilising"


###############################
##  Sampled trunc


sampled_trunc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled.WZA.csv")
hist(sampled_trunc_10000_s0.003[sampled_trunc_10000_s0.003$LA>0,]$LA, breaks = 20)
sampled_trunc_10000_s0.003_LA_thresh = 0.005
new_sampled_trunc_10000_s0.003 <- sampled_trunc_10000_s0.003 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_trunc_10000_s0.003_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled/"
sampled_trunc_10000_s0.003_suffix <- "_0.003_directional_d40n50_i10000.csv"
sampled_trunc_10000_s0.003_FP <- NumberOfTruePositives(new_sampled_trunc_10000_s0.003, sampled_trunc_10000_s0.003_path, sampled_trunc_10000_s0.003_suffix, sampled_trunc_10000_s0.003_LA_thresh)
sampled_trunc_10000_s0.003_FP$map <- "trunc Map"
sampled_trunc_10000_s0.003_FP$s_max  <- "0.003"
sampled_trunc_10000_s0.003_summary <- summarise( sampled_trunc_10000_s0.003_FP )
sampled_trunc_10000_s0.003_summary$map <- "Truncated"
sampled_trunc_10000_s0.003_summary$selection <- "Directional"

sampled_trunc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/trunc_sampled.WZA.csv")
hist(sampled_trunc_10000_Vs192[sampled_trunc_10000_Vs192$LA>0,]$LA, breaks = 20)
sampled_trunc_10000_Vs192_LA_thresh = 0.005
new_sampled_trunc_10000_Vs192 <- sampled_trunc_10000_Vs192 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(X_Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)
sampled_trunc_10000_Vs192_path <- "~/work/GEA/simulations/directionalSelection/G.2.3/Vs192/trunc_sampled/"
sampled_trunc_10000_Vs192_suffix <- "__0.5_192_d40n50_i10000.csv"
sampled_trunc_10000_Vs192_FP <- NumberOfTruePositives(new_sampled_trunc_10000_Vs192, sampled_trunc_10000_Vs192_path, sampled_trunc_10000_Vs192_suffix, sampled_trunc_10000_Vs192_LA_thresh)
sampled_trunc_10000_Vs192_FP$map <- "Truncated"
sampled_trunc_10000_Vs192_FP$s_max  <- "0.003"
sampled_trunc_10000_Vs192_summary <- summarise( sampled_trunc_10000_Vs192_FP )
sampled_trunc_10000_Vs192_summary$map <- "Truncated"
sampled_trunc_10000_Vs192_summary$selection <- "Stabilising"


stabilising_all <- rbind( bc_10000_Vs192_summary,
                  cline_10000_Vs192_summary,
                  trunc_10000_Vs192_summary)
stabilising_all$data <- "Complete"

stabilising_sampled <- rbind(sampled_bc_10000_Vs192_summary,
                             sampled_cline_10000_Vs192_summary,
                             sampled_trunc_10000_Vs192_summary)
stabilising_sampled$data <- "Sampled"

stabilising<- rbind( stabilising_all, stabilising_sampled)

directional_sampled <- rbind( sampled_bc_10000_s0.003_summary,
                              sampled_cline_10000_s0.003_summary,
                              sampled_trunc_10000_s0.003_summary)
directional_sampled$data <- "Sampled"
directional_all <- rbind( bc_10000_s0.003_summary,
                   cline_10000_s0.003_summary,
                   trunc_10000_s0.003_summary)
directional_all$data <- "Complete"
directional<- rbind( directional_all, directional_sampled)

#### Make a plot of the result

directionalPlot <- ggplot(data = directional, aes( x = top, y = mean, col = test, fill = test))+
  geom_line(lwd = 2)+
  geom_ribbon( aes( ymax = mean + se*2, ymin = mean - se*2, fill = test, colour=NULL), alpha = 0.3, lwd = 0)+
  theme_bw()+
  facet_grid(data~ map)+
  scale_colour_brewer('Test', palette = "Set2")+
  scale_fill_brewer('Test', palette = "Set2")+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )


stabilisingPlot <- ggplot(data = stabilising, aes( x = top, y = mean, col = test, fill = test))+
  geom_line(lwd = 2)+
  geom_ribbon( aes( ymax = mean + se*2, ymin = mean - se*2, fill = test, colour=NULL), alpha = 0.3, lwd = 0)+
  theme_bw()+
  facet_grid(data~ map)+
  scale_colour_brewer('Test', palette = "Set2")+
  scale_fill_brewer('Test', palette = "Set2")+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )


pdf("~/work/GEA/simulations/Plots/StabilisingSelection_dataComparison.pdf", width = 10, height = 6)
print(stabilisingPlot)
dev.off()

pdf("~/work/GEA/simulations/Plots/DirectionalSelection_dataComparison.pdf", width = 10, height = 6)
print(directionalPlot)
dev.off()



######## Plot effect size Distribution

bc_10000_Vs192$map <- "BC Map"
cline_10000_Vs192$map <- "Cline"
trunc_10000_Vs192$map <- "Truncated"
stablising_data <- rbind( bc_10000_Vs192, cline_10000_Vs192, trunc_10000_Vs192 )
stablising_data$model <- "Stabilising Selection"

bc_10000_s0.003$map <- "BC Map"
cline_10000_s0.003$map <- "Cline"
trunc_10000_s0.003$map <- "Truncated"
directional_data <- rbind( bc_10000_s0.003, cline_10000_s0.003, trunc_10000_s0.003 )
directional_data$model <- "Directional Selection"

effect_size_data <- rbind( stablising_data, directional_data)


effectSizeDistribution_dir <- ggplot()+
  geom_histogram( data=directional_data[directional_data$LA > 0, ], aes( x= LA ), bins = 50, fill = "#91cf60", col = "black")+
  geom_vline(xintercept = 0.005, col = "red", lty = 2)+
  facet_wrap(model~ map, scales = "free")+
  xlab("Covariance of Phenotype and Environment")+
  ylab("Count")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )

effectSizeDistribution_stab <- ggplot()+
  geom_histogram( data=stablising_data[stablising_data$LA > 0, ], aes( x= LA ), bins = 60, fill = "#fc8d59", col = "black")+
  geom_vline(xintercept = 0.005, col = "red", lty = 2)+
  facet_wrap(model~ map, scales = "free")+
    xlab("Covariance of Phenotype and Environment")+
  ylab("Count")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )


pdf("~/work/GEA/simulations/Plots/effectSizeDistributionPlot.pdf", width = 12, height = 5)
print(effectSizeDistribution)  
dev.off()
