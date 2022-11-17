rm(list = ls())

falsePositiveCounter <- function(gea_results){
  falsePositiveDF <- data.frame(0,0,0,0,0,0,0)
  names(falsePositiveDF) <-  c("top", "truePositive_WZA", "falsePositive_WZA",  "truePositive_WZA_BayPass", "falsePositive_WZA_BayPass", "truePositive_TC" , "falsePositive_TC")
  for (n in 1:200){
    
    gea_results_hits <- gea_results[gea_results$LA > 0.01,]
    gea_results_neutral <- gea_results[gea_results$distToPVE > 30000,]
    gea_results <- rbind(gea_results_hits, gea_results_neutral)
    
    WZA_sorted <- gea_results[order(gea_results$Z_kendall,decreasing = T),]
    WZA_BayPass_sorted <- gea_results[order(gea_results$Z_bayPass,decreasing = T),]
    top_candidate_sorted <- gea_results[order(gea_results$top_candidate_p,decreasing = F),]
    
    top_n_WZA <-WZA_sorted[1:n,]
    top_n_WZA_BayPass <-WZA_BayPass_sorted[1:n,]
    top_n_TC <-top_candidate_sorted[1:n,]
    truePositive_WZA = nrow( top_n_WZA[top_n_WZA$LA > 0.01,] )
    falsePositive_WZA = nrow( top_n_WZA) - truePositive_WZA 
    truePositive_WZA_BayPass = nrow( top_n_WZA_BayPass[top_n_WZA_BayPass$LA > 0.01,] )
    falsePositive_WZA_BayPass = nrow( top_n_WZA_BayPass) - truePositive_WZA_BayPass 
    truePositive_TC = nrow( top_n_TC[top_n_TC$LA > 0.01,] )
    falsePositive_TC = nrow( top_n_TC) - truePositive_TC 
    falsePositiveDF[n,] = c(n, truePositive_WZA, falsePositive_WZA,truePositive_WZA_BayPass, falsePositive_WZA_BayPass, truePositive_TC,falsePositive_TC)
  }
  return(falsePositiveDF)
}

falsePositiveCollater <- function(gea_dataFrame){
  gea_fp_list <- list()
  counter = 0
  for (Vs in c(100,200)){
    for (U_a in c(0.0001,0.001)){
      for (sig_a in c(0.3, 1)){
        print(c(Vs,U_a,sig_a))
        counter = counter +1
        gea_temp <- gea_dataFrame[( gea_dataFrame$Vs == Vs )&( gea_dataFrame$sig_a == sig_a )&( gea_dataFrame$U_a ==U_a  ),]
        fp_temp <- falsePositiveCounter(gea_temp)
        fp_temp$U_a <- U_a
        fp_temp$Vs <- Vs
        fp_temp$sig_a <- sig_a
        gea_fp_list[[counter]] <- fp_temp
        
      } 
    }
  }
  return (do.call("rbind",gea_fp_list))
  
}

labelIt <- function(FP_DF){
  FP_DF$sig_a <- factor(FP_DF$sig_a ,levels = c(0.3,1),
                  labels = c(expression(sigma[alpha]^2*" = 0.3" ),expression(sigma[alpha]^2*" = 1.0" )))

  FP_DF$U_a <- factor(FP_DF$U_a , levels = c(0.001,1e-04),
                labels = c(expression(U[alpha]*" = 0.001" ),expression(U[alpha]*" = 0.0001" )))


  FP_DF$Vs <- factor(FP_DF$Vs , levels = c(100,200),
               labels = c(expression(V[s]*" = 100" ),expression(V[s]*" = 200" )))
  return(FP_DF)
}

## READ IN ALL THE DATA
bc_map_all <- read.csv("/media/booker/HOWDY/GEA/G.1.1/FromCedar/BC_map_d40_n50_WZA_bayPass.csv")
bc_FP <- falsePositiveCollater(bc_map_all)

cline_map_all <- read.csv("/media/booker/HOWDY/GEA/G.1.1/FromCedar/cline_d40_n50_WZA_bayPass.csv")
cline_FP <- falsePositiveCollater(cline_map_all)

random_map_all <- read.csv("/media/booker/HOWDY/GEA/G.1.1/FromCedar/random_d40_n50_WZA_bayPass.csv")
random_FP <- falsePositiveCollater(random_map_all)

island_map_all <- read.csv("/media/booker/HOWDY/GEA/G.1.1/island_d40_n50_WZA.csv")
island_FP <- falsePositiveCollater(island_map_all)

library(ggplot2)

png("/media/booker/HOWDY/GEA/Plots/BC_Map_falsePositives.png", height = 12, width = 8, res = 200, units = "in")
ggplot( labelIt( bc_FP ), aes(x = top, y = truePositive_WZA))+
  geom_line()+
  geom_line(aes(x = top, y = truePositive_TC), col = "red")+
  geom_line(aes(x = top, y = truePositive_WZA_BayPass), col = "blue")+
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey")+
  scale_y_continuous("Number of True Positives",limits = c(0,200))+
  scale_x_continuous("Top n Hits",limits = c(0,200))+
  facet_grid( U_a + sig_a ~ Vs, labeller = label_parsed )+
  theme_bw()
dev.off()

png("/media/booker/HOWDY/GEA/Plots/cline_falsePositives.png", height = 12, width = 8, res = 200, units = "in")
ggplot(labelIt( cline_FP ), aes(x = top, y = truePositive_WZA))+
  geom_line()+
  geom_line(aes(x = top, y = truePositive_TC), col = "red")+
  geom_line(aes(x = top, y = truePositive_WZA_BayPass), col = "blue")+
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey")+
  scale_y_continuous("Number of True Positives",limits = c(0,200))+
  scale_x_continuous("Top n Hits",limits = c(0,200))+
  facet_grid(U_a+sig_a~Vs,  labeller = label_parsed)+
  theme_bw()
dev.off()

png("/media/booker/HOWDY/GEA/Plots/random_falsePositives.png", height = 12, width = 8, res = 200, units = "in")
ggplot(labelIt( random_FP ), aes(x = top, y = truePositive_WZA))+
  geom_line()+
  geom_line(aes(x = top, y = truePositive_TC), col = "red")+
  geom_line(aes(x = top, y = truePositive_WZA_BayPass), col = "blue")+
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "grey")+
  scale_y_continuous("Number of True",limits = c(0,200))+
  scale_x_continuous("Top n Hits",limits = c(0,200))+
  facet_grid(U_a+sig_a~Vs,  labeller = label_parsed)+
  theme_bw()
dev.off()

random_Vs100_Ua0.001_SigA0.3 <- random_map_all[(random_map_all$U_a == 0.001)&(random_map_all$Vs == 100)&(random_map_all$sig_a==0.3),]
random_Vs100_Ua0.001_SigA0.3$POS <- random_Vs100_Ua0.001_SigA0.3$position + 1e6*random_Vs100_Ua0.001_SigA0.3$rep

ggplot(data = random_Vs100_Ua0.001_SigA0.3[(random_Vs100_Ua0.001_SigA0.3$distToPVE > 25000)|(random_Vs100_Ua0.001_SigA0.3$LA > 0.01),], aes( x = POS, y = Z_bayPass, col = LA> 0))+
  geom_point()
ggplot(data = random_Vs100_Ua0.001_SigA0.3[(random_Vs100_Ua0.001_SigA0.3$distToPVE > 25000)|(random_Vs100_Ua0.001_SigA0.3$LA > 0.01),], aes( x = POS, y = Z_kendall, col = LA> 0))+
  geom_point()

