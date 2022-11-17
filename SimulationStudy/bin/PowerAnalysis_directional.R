rm(list = ls()) 

library(ggplot2)
library(reshape2)
library(dplyr)

pal <- c("#d7191c","#fdae61","#abd9e9","#2c7bb6")

powerFuncParametric <- function(gea_data){
  
  gea_data <-   gea_data[  -log10(gea_data$LA) > 100,]

  gea_data$TC_empirical_p <-  rank( gea_data$top_candidate_p) /length(gea_data$top_candidate_p) 

  #  gea_data$top_candidate_hist <- gea_data$top_candidate_p < 0.001
  gea_data$Z_hits <- gea_data$Z > qnorm(1 - 0.01, mean =mean(gea_data$Z),sd =  sd(gea_data$Z))
  gea_data$top_candidate_hist <- gea_data$TC_empirical_p < 0.01
  #gea_data$top_candidate_hist <- gea_data$TC_empirical_p < sum(gea_data$Z_hits)/length(gea_data$TC_empirical_p)
  
  gea_data$big_hitter <- gea_data$bigHits >0
  #  gea_data$Z_hits <- gea_data$Z > quantile(gea_data$Z, 0.99)
  #  gea_data$Z_hits<-gea_data$Z%in%gea_data$Z[order(gea_data$Z, decreasing = T)][1:sum(gea_data$top_candidate_p < 0.001)]
  
  
  mod1 <- glm(top_candidate_hist ~ LA, family = binomial, data = gea_data)
  mod2 <- glm(Z_hits ~ LA, family = binomial, data = gea_data)
  mod3 <- glm(big_hitter ~ LA, family = binomial, data = gea_data)
  
  x<-1:1000/1000
  
  y1<-predict(mod1,  list(LA = x), type = "response")
  y2<-predict(mod2,  list(LA = x), type = "response")
  y3<-predict(mod3,  list(LA = x), type = "response")
  
  return(data.frame(PVE = x, TC=y1, WZA=y2, snps = y3))
}

powerFuncNonParametric <- function(gea_data, p_value){
## Let's look at the probabilityies that genes are in the top XXX of the data
  
  gea_data$TC_empirical_p <-  rank( gea_data$top_candidate_p) /length(gea_data$top_candidate_p) 
  gea_data$Z_empirical_p <-  (length(gea_data$Z) - rank( gea_data$Z)) /length(gea_data$Z) 
  
  gea_data$LA_index <- -log10(gea_data$LA)
  gea_data <-   gea_data[  gea_data$LA_index > 3,]
  
  
  p_value <- p_value
  
  gea_data$top_candidate_hits <- gea_data$TC_empirical_p < p_value
  gea_data$Z_hits <- gea_data$Z_empirical_p < p_value
  
  mod1 <- glm(top_candidate_hits ~ LA_index, family = binomial, data = gea_data)
  mod2 <- glm(Z_hits ~ LA_index, family = binomial, data = gea_data)
  
  
  x<-1:3000/100
  
  y1<-predict(mod1,  list(LA_index = x), type = "response")
  y2<-predict(mod2,  list(LA_index = x), type = "response")
  
  return(data.frame(PVE = x, TC=y1, WZA=y2))
}



correctTC <- function(df){
  q <- df[df$hits!=0,]
  
  expected <- sum(q$hits)/sum(q$SNPs)
  
  return(dbinom(df$hits, df$SNPs, expected))
  
}


powerFuncNonParametricPerRep <- function(gea_data, p_value, LA_threshold){
  ## Let's look at the probabilityies that genes are in the top XXX of the data
  
  #  gea_data$TC_empirical_p <-  rank( gea_data$top_candidate_p) /length(gea_data$top_candidate_p) 
  #  gea_data$Z_empirical_p <-  (length(gea_data$Z) - rank( gea_data$Z)) /length(gea_data$Z) 

#  gea_data <- bc  
  gea_data$LA_index <- gea_data$LA
  gea_data <-   gea_data[  gea_data$LA_index > LA_threshold,]
  
  gea_data$top_candidate_hits <- gea_data$top_candidate_p < p_value
  gea_data$Z_hits <- gea_data$Z_empirical_p < p_value
  
  mod1 <- glm(top_candidate_hits ~ LA_index, family = binomial, data = gea_data)
  mod2 <- glm(Z_hits ~ LA_index, family = binomial, data = gea_data)
  
  
  x<-1:3000/100
  
  y1<-predict(mod1,  list(LA_index = x), type = "response")
  y2<-predict(mod2,  list(LA_index = x), type = "response")
  
  return(data.frame(PVE = x, TC=y1, WZA=y2))
}


falsePositiveCounter <- function(gea_results, LA_threshold){
  falsePositiveDF <- data.frame(0,0,0,0)
  names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC")
  for (rep in 0:19){
    temp <- gea_results[gea_results$rep == rep,]
    numberLA <- length(unname(temp[temp$LA > LA_threshold,]$gene))
    for (n in 0:49){

      WZA_slice <- temp[temp$Z_empirical_p <= (n+1)/1000,]
      TC_slice <- temp[temp$TC_empirical_p <= (n+1)/1000,]
      
      truePositive_WZA = sum( -log10(WZA_slice$LA) > LA_threshold )/numberLA
      truePositive_TC = sum( -log10(TC_slice$LA) > LA_threshold )/numberLA
      falsePositiveDF[(n*50)+rep,] <- c(rep+1, n+1, truePositive_WZA, truePositive_TC)
    }
  }
  return(falsePositiveDF)
}

############################################
############################################
############################################
###
###
### Sampled data - Power
###
###
############################################
############################################
############################################

bc_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_10_sampled.csv")
bc_10$top_candidate_p<- correctTC(bc_10)

str(bc_10)
temp <- bc_10[bc_10$rep == 1,]


newBC_10 <- bc_10 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_10_results_0.05 <- powerFuncNonParametricPerRep(newBC_10, 0.05)
bc_10_results_0.05$alpha <- 0.05
bc_10_results_0.01 <- powerFuncNonParametricPerRep(newBC_10, 0.01)
bc_10_results_0.01$alpha <- 0.01
bc_10_results <- rbind(bc_10_results_0.01,bc_10_results_0.05)
bc_10_results$alpha <-as.factor(bc_10_results$alpha)

bc_10_meltedPowerSet<- melt(bc_10_results, id = c("PVE", "alpha"))
bc_10_meltedPowerSet$variable <- factor(bc_10_meltedPowerSet$variable , 
                                  levels = c("TC","WZA"),
                                  labels = c("Top-Candidate Test", "The WZA"))

bc_10_meltedPowerSet$alpha <- factor(bc_10_meltedPowerSet$alpha , 
                                  levels = c(0.01,0.05),
                                  labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

bc_100 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_100_sampled.csv"))
bc_100$top_candidate_p<- correctTC(bc_100)

newBC_100 <- bc_100 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_100_results_0.05 <- powerFuncNonParametricPerRep(newBC_100, 0.05)
bc_100_results_0.05$alpha <- 0.05
bc_100_results_0.01 <- powerFuncNonParametricPerRep(newBC_100, 0.01)
bc_100_results_0.01$alpha <- 0.01
bc_100_results <- rbind(bc_100_results_0.01,bc_100_results_0.05)
bc_100_results$alpha <-as.factor(bc_100_results$alpha)

bc_100_meltedPowerSet<- melt(bc_100_results, id = c("PVE", "alpha"))
bc_100_meltedPowerSet$variable <- factor(bc_100_meltedPowerSet$variable , 
                                        levels = c("TC","WZA"),
                                        labels = c("Top-Candidate Test", "The WZA"))

bc_100_meltedPowerSet$alpha <- f}actor(bc_100_meltedPowerSet$alpha , 
                                     levels = c(0.01,0.05),
                                     labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

bc_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_1000_sampled.csv")
bc_1000$top_candidate_p<- correctTC(bc_1000)

newBC_1000 <- bc_1000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_1000_results_0.05 <- powerFuncNonParametricPerRep(newBC_1000, 0.05)
bc_1000_results_0.05$alpha <- 0.05
bc_1000_results_0.01 <- powerFuncNonParametricPerRep(newBC_1000, 0.01)
bc_1000_results_0.01$alpha <- 0.01
bc_1000_results <- rbind(bc_1000_results_0.01,bc_1000_results_0.05)
bc_1000_results$alpha <-as.factor(bc_1000_results$alpha)

bc_1000_meltedPowerSet<- melt(bc_1000_results, id = c("PVE", "alpha"))
bc_1000_meltedPowerSet$variable <- factor(bc_1000_meltedPowerSet$variable , 
                                        levels = c("TC","WZA"),
                                        labels = c("Top-Candidate Test", "The WZA"))

bc_1000_meltedPowerSet$alpha <- factor(bc_1000_meltedPowerSet$alpha , 
                                     levels = c(0.01,0.05),
                                     labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

  bc_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map.WZA.sampled.csv")
bc_10000$top_candidate_p<- correctTC(bc_10000)

hist(bc_10000[bc_10000$LA!=0,]$LA)

newBC_10000 <- bc_10000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_10000_results_0.05 <- powerFuncNonParametricPerRep(newBC_10000, 0.05, 1)
bc_10000_results_0.05$alpha <- 0.05
bc_10000_results_0.01 <- powerFuncNonParametricPerRep(newBC_10000, 0.01, 1)
bc_10000_results_0.01$alpha <- 0.01
bc_10000_results <- rbind(bc_10000_results_0.01,bc_10000_results_0.05)
bc_10000_results$alpha <-as.factor(bc_10000_results$alpha)

bc_10000_meltedPowerSet<- melt(bc_10000_results, id = c("PVE", "alpha"))
bc_10000_meltedPowerSet$variable <- factor(bc_10000_meltedPowerSet$variable , 
                                        levels = c("TC","WZA"),
                                        labels = c("Top-Candidate Test", "The WZA"))

bc_10000_meltedPowerSet$alpha <- factor(bc_10000_meltedPowerSet$alpha , 
                                     levels = c(0.01,0.05),
                                     labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))


############################################
############################################


bc_10_meltedPowerSet$blockSize <- 10
bc_100_meltedPowerSet$blockSize <- 100
bc_1000_meltedPowerSet$blockSize <- 1000
bc_10000_meltedPowerSet$blockSize <- 10000

meltedResults <- rbind( bc_10_meltedPowerSet, bc_100_meltedPowerSet, bc_1000_meltedPowerSet, bc_10000_meltedPowerSet)
meltedResults <- rbind(bc_10000_meltedPowerSet)

ggplot(data = meltedResults, aes(y = value, x = PVE, col = variable))+
  geom_line(lwd = 1)+
  facet_grid(alpha~blockSize, labeller = label_parsed)+
  scale_colour_brewer(palette = "Set2")+
  scale_y_continuous("Power")+
  scale_x_continuous("LA Index")+
  theme_bw()+
  theme(
    legend.title = element_blank()
)



############################################
############################################
############################################
###
###
### Sampled data - False Positives
###
###
############################################
############################################
############################################

bc_10_FP <- na.omit( falsePositiveCounter(newBC_10, 1) )
bc_10_FP$blockSize <- 10
bc_1000_FP <- na.omit( falsePositiveCounter(newBC_1000, 1) )
bc_1000_FP$blockSize <- 1000
bc_10000_FP <- na.omit( falsePositiveCounter(newBC_10000, 1) )
bc_10000_FP$blockSize <- 10000

FP <- rbind(bc_10_FP,bc_1000_FP,bc_10000_FP)
FP <- rbind(bc_10000_FP)

ggplot(data = FP, aes( x = top, y = truePositive_WZA))+
  geom_smooth(aes(col = "WZA", fill = "WZA"))+
  geom_smooth(aes( x = top, y = truePositive_TC, col = "Top-Candidate Test", fill = "Top-Candidate Test"))+
  theme_bw()+
  facet_grid(blockSize ~.)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # genes")



############################################
############################################
############################################
###
###
### Unsampled data - Power
###
###
############################################
############################################
############################################



bc_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_10.csv")
bc_10$top_candidate_p<- correctTC(bc_10)

newBC_10 <- bc_10 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_10_results_0.05 <- powerFuncNonParametricPerRep(newBC_10, 0.05)
bc_10_results_0.05$alpha <- 0.05
bc_10_results_0.01 <- powerFuncNonParametricPerRep(newBC_10, 0.01)
bc_10_results_0.01$alpha <- 0.01
bc_10_results <- rbind(bc_10_results_0.01,bc_10_results_0.05)
bc_10_results$alpha <-as.factor(bc_10_results$alpha)

bc_10_meltedPowerSet<- melt(bc_10_results, id = c("PVE", "alpha"))
bc_10_meltedPowerSet$variable <- factor(bc_10_meltedPowerSet$variable , 
                                        levels = c("TC","WZA"),
                                        labels = c("Top-Candidate Test", "The WZA"))

bc_10_meltedPowerSet$alpha <- factor(bc_10_meltedPowerSet$alpha , 
                                     levels = c(0.01,0.05),
                                     labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

bc_100 <- na.omit(read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_100.csv"))
bc_100$top_candidate_p<- correctTC(bc_100)

newBC_100 <- bc_100 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_100_results_0.05 <- powerFuncNonParametricPerRep(newBC_100, 0.05)
bc_100_results_0.05$alpha <- 0.05
bc_100_results_0.01 <- powerFuncNonParametricPerRep(newBC_100, 0.01)
bc_100_results_0.01$alpha <- 0.01
bc_100_results <- rbind(bc_100_results_0.01,bc_100_results_0.05)
bc_100_results$alpha <-as.factor(bc_100_results$alpha)

bc_100_meltedPowerSet<- melt(bc_100_results, id = c("PVE", "alpha"))
bc_100_meltedPowerSet$variable <- factor(bc_100_meltedPowerSet$variable , 
                                         levels = c("TC","WZA"),
                                         labels = c("Top-Candidate Test", "The WZA"))

bc_100_meltedPowerSet$alpha <- factor(bc_100_meltedPowerSet$alpha , 
                                      levels = c(0.01,0.05),
                                      labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

bc_1000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_1000.csv")
bc_1000$top_candidate_p<- correctTC(bc_1000)

newBC_1000 <- bc_1000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_1000_results_0.05 <- powerFuncNonParametricPerRep(newBC_1000, 0.05)
bc_1000_results_0.05$alpha <- 0.05
bc_1000_results_0.01 <- powerFuncNonParametricPerRep(newBC_1000, 0.01)
bc_1000_results_0.01$alpha <- 0.01
bc_1000_results <- rbind(bc_1000_results_0.01,bc_1000_results_0.05)
bc_1000_results$alpha <-as.factor(bc_1000_results$alpha)

bc_1000_meltedPowerSet<- melt(bc_1000_results, id = c("PVE", "alpha"))
bc_1000_meltedPowerSet$variable <- factor(bc_1000_meltedPowerSet$variable , 
                                          levels = c("TC","WZA"),
                                          labels = c("Top-Candidate Test", "The WZA"))

bc_1000_meltedPowerSet$alpha <- factor(bc_1000_meltedPowerSet$alpha , 
                                       levels = c(0.01,0.05),
                                       labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))

############################################
############################################

bc_10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_10000.csv")
bc_10000$top_candidate_p<- correctTC(bc_10000)

newBC_10000 <- bc_10000 %>%
  group_by(rep) %>%
  mutate(Z_empirical_p = order(order(Z, decreasing = T))/ 1000,
         TC_empirical_p = order(order(-log(top_candidate_p), decreasing = T))/ 1000)


bc_10000_results_0.05 <- powerFuncNonParametricPerRep(newBC_10000, 0.05)
bc_10000_results_0.05$alpha <- 0.05
bc_10000_results_0.01 <- powerFuncNonParametricPerRep(newBC_10000, 0.01)
bc_10000_results_0.01$alpha <- 0.01
bc_10000_results <- rbind(bc_10000_results_0.01,bc_10000_results_0.05)
bc_10000_results$alpha <-as.factor(bc_10000_results$alpha)

bc_10000_meltedPowerSet<- melt(bc_10000_results, id = c("PVE", "alpha"))
bc_10000_meltedPowerSet$variable <- factor(bc_10000_meltedPowerSet$variable , 
                                           levels = c("TC","WZA"),
                                           labels = c("Top-Candidate Test", "The WZA"))

bc_10000_meltedPowerSet$alpha <- factor(bc_10000_meltedPowerSet$alpha , 
                                        levels = c(0.01,0.05),
                                        labels = c(expression(alpha*" = 0.01"), expression(alpha*" = 0.05")))


############################################
############################################


bc_10_meltedPowerSet$blockSize <- 10
#bc_100_meltedPowerSet$blockSize <- 100
bc_1000_meltedPowerSet$blockSize <- 1000
bc_10000_meltedPowerSet$blockSize <- 10000

meltedResults <- rbind( bc_10_meltedPowerSet, bc_1000_meltedPowerSet, bc_10000_meltedPowerSet)

ggplot(data = meltedResults, aes(y = value, x = PVE, col = variable))+
  geom_line(lwd = 1)+
  facet_grid(alpha~blockSize, labeller = label_parsed)+
  scale_colour_brewer(palette = "Set2")+
  scale_y_continuous("Power")+
  scale_x_continuous("-log10(p-value)")+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )


############################################
############################################
############################################
###
###
### Sample Data False Positives
###
###
############################################
############################################
############################################


bc_10_FP <- na.omit( falsePositiveCounter(newBC_10) )
bc_10_FP$blockSize <- 10
bc_1000_FP <- na.omit( falsePositiveCounter(newBC_1000) )
bc_1000_FP$blockSize <- 1000
bc_10000_FP <- na.omit( falsePositiveCounter(newBC_10000, 1) )
bc_10000_FP$blockSize <- 10000

FP <- rbind(bc_10_FP,bc_1000_FP,bc_10000_FP)
FP <- bc_10000_FP
ggplot(data = FP, aes( x = top, y = truePositive_WZA))+
  geom_smooth(aes(col = "WZA", fill = "WZA"))+
  geom_smooth(aes( x = top, y = truePositive_TC, col = "Top-Candidate Test", fill = "Top-Candidate Test"))+
  theme_bw()+
  facet_grid(blockSize ~.)+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # genes")


bc_10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BC_Map_recombinationVariation/BC_Map_0.003_10_sampled/1_0.003_1.12Loci.directionalSelection_10bp_p40_n50.csv")


gea_results <- newBC_10000


path <- "~/work/GEA/simulations/directionalSelection/G.2.3/s0.0055/BC_Map_sampled//"
suffix <- "_0.0055_1.12Loci.directionalSelection_d40_n50.csv"
falsePositiveDF <- data.frame(0,0,0,0,0)
names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC","truePositive_SNPs")

for (rep in 1:20){

  temp <- gea_results[gea_results$rep == rep,]
  numberLA <- length(droplevels(temp[temp$LA >1,]$gene))
  LA_genes <- droplevels(unique(temp[temp$LA > 1,]$gene))
  
  snps <- read.csv(paste(path, rep, suffix, sep = ''))

  snps <- snps[ with(snps, order(ave(pop_k_tau_p_value, gene, FUN = min), pop_k_tau_p_value)),]
  snp_based_genes <- unique( snps$gene )
   
  for (n in 1:50){

      WZA_slice <- temp[temp$Z_empirical_p <= (n+1)/1000,]
      TC_slice <- temp[temp$TC_empirical_p <= (n+1)/1000,]
      SNP_slice <- snp_based_genes[1:n]
      
#      SNP_slice <- snps[snps$rank <= n,] ## Grab the top n SNPs from the genome scan
            
#      numberLA_SNP_slice <- length(unname(unique(SNP_slice[-log10(SNP_slice$LA) > 3.0,]$gene)))

      truePositive_WZA = sum( WZA_slice$LA > 1 )/numberLA
      truePositive_TC = sum( TC_slice$LA > 1 )/numberLA
      truePositive_SNPs = sum(SNP_slice%in%LA_genes)/numberLA
      
      falsePositiveDF[(50*(rep-1))+n,] <- c(rep, n, truePositive_WZA, truePositive_TC, truePositive_SNPs)

    }
}


ggplot(data = falsePositiveDF, aes( x = top, y = truePositive_WZA))+
  geom_jitter(size = 0, alpha = 0)+
  stat_summary(data = falsePositiveDF, aes(y = truePositive_WZA, x = top, group = 1, col = "WZA"), fun.y = mean, geom = "line", group = 1, lwd = 2)+
  stat_summary(data = falsePositiveDF, aes(y = truePositive_TC, x = top, group = 2 , col = "Top-Candidate Test"), fun.y = mean, geom = "line", group = 2, lwd = 2)+
  stat_summary(data = falsePositiveDF, aes(y = truePositive_SNPs, x = top, group = 3, col = "SNP-Based Test"), fun.y = mean, geom = "line", group = 3, lwd = 2)+
  theme_bw()+
  scale_colour_brewer('', palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = F)+
  scale_y_continuous("Proportion of True Positives Detected")+
  scale_x_continuous("Top # Genes")



falsePositiveDF[  falsePositiveDF$rep == 15,]

XXX<-snps[ sample(  10000, 1000 ),]
XXX$gene <- droplevels(XXX$gene)
df <- XXX[ with(XXX, order(ave(pop_k_tau_p_value, gene, FUN = min), pop_k_tau_p_value)),]

unique( df$gene )
