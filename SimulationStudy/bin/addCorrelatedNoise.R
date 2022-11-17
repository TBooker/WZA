rm(list = ls())


library(ggplot2)
library(ggpubr)
library(reshape2)
## Sample correlated values for the population environments

# Here's a function to simulate correlated values'

getBiCop <- function(x, rho){
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  C <- chol(C)

  X2 <- rnorm(length(x))
  X <- cbind(x, X2)
  
  df <- X %*% C
  
  return(df)
      
}





## Now, read in the BC map environments

BC_envVec <- read.csv("~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt", header = F)$V1
## Set the seed so we alsways get the same result for these
set.seed(10)
BC_envVec_rho0.1 <- data.frame(getBiCop( BC_envVec, 0.1))
write.table( BC_envVec_rho0.1$X2, "~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.cor0.1.txt", row.names = F, col.names = F)
BC_envVec_rho0.1$correlation <- "Correlation = 0.1"

set.seed(10)
BC_envVec_rho0.3 <- data.frame(getBiCop( BC_envVec, 0.3))
write.table( BC_envVec_rho0.3$X2, "~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.cor0.3.txt", row.names = F, col.names = F)
BC_envVec_rho0.3$correlation <- "Correlation = 0.3"

set.seed(10)
BC_envVec_rho0.5 <- data.frame(getBiCop( BC_envVec, 0.5))
write.table( BC_envVec_rho0.5$X2, "~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.cor0.5.txt", row.names = F, col.names = F)
BC_envVec_rho0.5$correlation <- "Correlation = 0.5"

set.seed(10)
BC_envVec_rho0.8 <- data.frame(getBiCop( BC_envVec, 0.8))
write.table( BC_envVec_rho0.8$X2, "~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.cor0.8.txt", row.names = F, col.names = F)
BC_envVec_rho0.8$correlation <- "Correlation = 0.8"

correlated_environments <- rbind(BC_envVec_rho0.1,
                                 BC_envVec_rho0.3,
                                 BC_envVec_rho0.5,
                                 BC_envVec_rho0.8)



stab_bc_map_cor_1.0 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.summary.csv")
stab_bc_map_cor_1.0$cor <- "Correlation = 1.0"

stab_bc_map_cor_0.8 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled_cor0.8.WZA.summary.csv")
stab_bc_map_cor_0.8$cor <- "Correlation = 0.8"

stab_bc_map_cor_0.5 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled_cor0.5.WZA.summary.csv")
stab_bc_map_cor_0.5$cor <- "Correlation = 0.5"

stab_bc_map_cor_0.3 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled_cor0.3.WZA.summary.csv")
stab_bc_map_cor_0.3$cor <- "Correlation = 0.3"

stab_bc_map_cor_0.1 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled_cor0.1.WZA.summary.csv")
stab_bc_map_cor_0.1$cor <- "Correlation = 0.1"

Stabilising_cors<- rbind(stab_bc_map_cor_0.1,
                         stab_bc_map_cor_0.3,
                         stab_bc_map_cor_0.5,
                         stab_bc_map_cor_0.8,
                         stab_bc_map_cor_1.0)

Stabilising_cors<- rbind(stab_bc_map_cor_0.1,
                         stab_bc_map_cor_0.3,
                         stab_bc_map_cor_0.5,
                         stab_bc_map_cor_0.8)



Stabilising_id_vars <- names(Stabilising_cors)[!(names(Stabilising_cors)%in%c( "PCE_SNP_uncor", "PCE_WZA_empR", "PCE_SNP_bayP", "PCE_WZA_bayP"))]
Stabilising_TruePositives <- melt( Stabilising_cors, id = Stabilising_id_vars  )
Stabilising_TruePositives$data <- "Uncorrected p-values"
Stabilising_TruePositives$variable <- factor(Stabilising_TruePositives$variable,
                                                   levels = c( "PCE_WZA_empR",  "PCE_SNP_uncor", "PCE_WZA_bayP", "PCE_SNP_bayP"))





paletteCB <- c("#D55E00",
               "#999999",
               "#CC79A7",
               "#009E73")
library(scales)

labs1 = c(expression("WZA"*tau),
          expression("Kendall's "*tau),
          expression("WZA"[BP]),
          expression("BayPass"))



bc_plot<- ggplot()+
  geom_line(data = Stabilising_TruePositives[Stabilising_TruePositives$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 0.9)+
  facet_wrap(. ~ cor)+
  theme_bw()+
  scale_color_manual("Method",
                     values = paletteCB,
                     labels = labs1)+  scale_y_continuous(expression("Proportion of Cov["*italic("phen, env")*"] explained"), expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 30, size = 8),
    strip.text.x =  element_text(size = 10),
    legend.text.align = 0,   
    plot.title = element_text(hjust = 0.5, size = 15)
  )



library(ggplot2)
library(ggpubr)

env_plot <- ggplot(data = correlated_environments , aes( x= X1, y = X2))+
  geom_point(pch = 21, fill = "grey")+
  theme_bw()+
  facet_grid(~correlation)+
  scale_y_continuous("Measured Environment (Arbitrary Units)")+
  scale_x_continuous("Phenotypic Optimum (Aribitrary Units)")+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


pdf("~/work/GEA/simulations/Plots/correlatedEnvironments.pdf", height = 4, width = 7)
print( env_plot )
dev.off()

pdf("~/work/GEA/simulations/Plots/correlatedEnvironments_BCmapResults.pdf", height = 4, width = 4.5)
print( bc_plot )
dev.off()




CVE_Directional_cors_TruePositives_Uncor <-  melt( Directional_cors, 
                                               id = c( "X" , "top" ,"FD_SNP_bayP", 
                                                       "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor",
                                                       "FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR", 
                                                       "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_uncor",
                                                       "TP_WZA_empR","TP_SNP_bayP", "TP_TC_bayP",
                                                       "TP_WZA_bayP" ,"PCE_WZA_uncor", "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP",
                                                       "rep", "cor")  )
CVE_Directional_cors_TruePositives_Uncor$data <- "Uncorrected GEA"
CVE_Directional_cors_TruePositives_Uncor$variable <- factor(CVE_Directional_cors_TruePositives_Uncor$variable,
                                                        levels = c( "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR"),
                                                        labels = c( "SNP-based", "Top-Candidate", "WZA"))


CVE_Directional_cors_TruePositives_Baypass <- melt( Directional_cors, 
                                                    id = c( "X" , "top" ,"FD_SNP_bayP", 
                                                            "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor",
                                                            "FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR", 
                                                            "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_uncor",
                                                            "TP_WZA_empR","TP_SNP_bayP", "TP_TC_bayP",
                                                            "TP_WZA_bayP" ,"PCE_WZA_uncor",
                                                            "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR",
                                                            "rep", "cor")  )
                                                    
                                                    
CVE_Directional_cors_TruePositives_Baypass$data <-  "BayPass"
CVE_Directional_cors_TruePositives_Baypass$variable <- factor(CVE_Directional_cors_TruePositives_Baypass$variable,
                                                          levels = c( "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP"),
                                                          labels = c( "SNP-based", "Top-Candidate", "WZA"))



bc_plot<- ggplot()+
  geom_line(data = CVE_Directional_cors_TruePositives_Uncor[CVE_Directional_cors_TruePositives_Uncor$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 0.9)+
  geom_line(data = CVE_Directional_cors_TruePositives_Baypass[CVE_Directional_cors_TruePositives_Baypass$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ cor)+
  theme_bw()+
  scale_color_manual("Method", values = paletteCB)+
  scale_y_continuous("Proportion of Cov(Phen.,Env.) explained", expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 30),
    strip.text.x =  element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


library(ggplot2)
library(ggpubr)

env_plot <- ggplot(data = correlated_environments , aes( x= X1, y = X2))+
  geom_point(pch = 21, fill = "grey")+
  theme_bw()+
  facet_grid(~correlation)+
  scale_y_continuous("Measured Environment (Arbitrary Units)")+
  scale_x_continuous("Phenotypic Optimum (Aribitrary Units)")+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


pdf("~/work/GEA/simulations/Plots/correlatedEnvironments.pdf", height = 4, width = 7)
print( env_plot )
dev.off()

#pdf("~/work/GEA/simulations/Plots/correlatedEnvironments_BCmapResults.pdf", height = 4, width = 8.5)
#print( bc_plot )
#dev.off()



Directional_cors<- rbind(stab_bc_map_cor_0.5,
                         stab_bc_map_cor_0.8,
                         stab_bc_map_cor_1.0)

Directional_cors_TruePositives_Uncor_ID_vars <- names(Directional_cors)[!(names(Directional_cors)%in%c( "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR"))]
Directional_cors_TruePositives_Uncor <-  melt( Directional_cors, id = Directional_cors_TruePositives_Uncor_ID_vars  )
Directional_cors_TruePositives_Uncor$data <- "Uncorrected p-values"
Directional_cors_TruePositives_Uncor$variable <- factor(Directional_cors_TruePositives_Uncor$variable,
                                                        levels = c( "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR"),
                                                        labels = c( "SNP-based", "Top-Candidate", "WZA"))


Directional_cors_TruePositives_Baypass_ID_vars <- names(Directional_cors)[!(names(Directional_cors)%in% c( "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP"))]
Directional_cors_TruePositives_Baypass <- melt( Directional_cors, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor",  "rep", "cor")  )
Directional_cors_TruePositives_Baypass$data <-  "BayPass"
Directional_cors_TruePositives_Baypass$variable <- factor(Directional_cors_TruePositives_Baypass$variable,
                                                          levels = c( "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP"),
                                                          labels = c( "SNP-based", "Top-Candidate", "WZA"))


talk_plot<- ggplot()+
  geom_line(data = Directional_cors_TruePositives_Uncor[Directional_cors_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_cors_TruePositives_Baypass[Directional_cors_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ cor)+
  theme_bw()+
  scale_color_manual("Method", values = paletteCB)+
  scale_y_continuous(expression("Proportion of Cov["*italic("phen, env")*"] explained"), expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 30),
    strip.text.x =  element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 15)
  )



pdf("~/work/GEA/simulations/Plots/correlatedEnvironments_BCmapResults_forTalk.pdf", height = 6, width = 7)
print( talk_plot )
dev.off()
