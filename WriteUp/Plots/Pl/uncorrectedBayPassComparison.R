rm(list = ls())

## Read in the necessary packages
library(reshape2)
library(ggplot2)
library(ggpubr)

dir_bc_map <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.summary.csv")
dir_bc_map$map <- "BC Map"
  
dir_cline <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/cline_sampled.WZA.summary.csv")
dir_cline$map <- "Gradient"

dir_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.summary.csv")
dir_trunc$map <- "Truncated"

Directional<- rbind(dir_bc_map, dir_cline, dir_trunc)

Directional_TruePositives_Uncor <- melt( Directional, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP" ,"rep", "map")  )
Directional_TruePositives_Uncor$data <- "Uncorrected p-values"
Directional_TruePositives_Uncor$variable <- factor(Directional_TruePositives_Uncor$variable,
                                                   levels = c( "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_uncor", "TP_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA - p-values", "WZA - empRho"))

Directional_TruePositives_Baypass <- melt( Directional, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR", "rep", "map")  )
Directional_TruePositives_Baypass$data <-  "BayPass"
Directional_TruePositives_Baypass$variable <- factor(Directional_TruePositives_Baypass$variable,
                                                     levels = c( "TP_SNP_bayP", "TP_TC_bayP", "TP_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA - p-values"))



stab_bc_map <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.summary.csv")
stab_bc_map$map <- "BC Map"

stab_cline <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/cline_sampled.WZA.summary.csv")
stab_cline$map <- "Gradient"

stab_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/trunc_sampled.WZA.summary.csv")
stab_trunc$map <- "Truncated"

Stabilising<- rbind(stab_bc_map, stab_cline, stab_trunc)

Stabilising_TruePositives_Uncor <- melt( Stabilising, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP" ,"rep", "map")  )
Stabilising_TruePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_TruePositives_Uncor$variable <- factor(Stabilising_TruePositives_Uncor$variable,
                                                   levels = c( "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_uncor", "TP_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA - p-values", "WZA - empRho"))

Stabilising_TruePositives_Baypass <- melt( Stabilising, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR", "rep", "map")  )
Stabilising_TruePositives_Baypass$data <-  "BayPass"
Stabilising_TruePositives_Baypass$variable <- factor(Stabilising_TruePositives_Baypass$variable,
                                                     levels = c( "TP_SNP_bayP", "TP_TC_bayP", "TP_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA - p-values"))





Directional_TruePositives_Uncor$selection <- "Directional Selection"
Stabilising_TruePositives_Uncor$selection <- "Stabilising Selection"

TruePositives_Uncor <- rbind(Directional_TruePositives_Uncor, Stabilising_TruePositives_Uncor)

Directional_TruePositives_Baypass$selection <- "Directional Selection"
Stabilising_TruePositives_Baypass$selection <- "Stabilising Selection"

TruePositives_Baypass <- rbind(Directional_TruePositives_Baypass, Stabilising_TruePositives_Baypass)

paletteCB <- c("#ca0020",
               "#f4a582",
               "#92c5de",
               "#0571b0")

directional_truePositives <- ggplot()+
  geom_line(data = Directional_TruePositives_Uncor[Directional_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_TruePositives_Baypass[Directional_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Directional Selection")+
  scale_color_manual("Test", values = paletteCB)+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


stablising_truePositives <- ggplot()+
  geom_line(data = Stabilising_TruePositives_Uncor[Stabilising_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Stabilising_TruePositives_Baypass[Stabilising_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Stablising Selection")+
  scale_color_manual("Test", values = paletteCB)+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
    
  )

combinedPlot <- ggarrange(directional_truePositives,stablising_truePositives, ncol= 2, common.legend = T, legend = "right", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_TruePositives.pdf", height = 9.5, width = 7)
print(combinedPlot)
dev.off()



### Now do the false positives



Directional_FalsePositives_Uncor <- melt( Directional, id = c( "X" , "top" ,"FD_SNP_bayP", "TP_SNP_uncor",  "FD_TC_bayP", "TP_TC_uncor","FD_WZA_bayP","TP_WZA_uncor", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP" ,"rep", "map")  )
Directional_FalsePositives_Uncor$data <- "Uncorrected p-values"
Directional_FalsePositives_Uncor$variable <- factor(Directional_FalsePositives_Uncor$variable,
                                                   levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_uncor"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA"))

Directional_FalsePositives_Baypass <- melt( Directional, id = c( "X" , "top" ,"TP_SNP_bayP", "FD_SNP_uncor",  "TP_TC_bayP", "FD_TC_uncor","TP_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor",  "rep", "map")  )
Directional_FalsePositives_Baypass$data <-  "BayPass"
Directional_FalsePositives_Baypass$variable <- factor(Directional_FalsePositives_Baypass$variable,
                                                     levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA"))




Stabilising_FalsePositives_Uncor <- melt( Stabilising, id = c( "X" , "top" ,"FD_SNP_bayP", "TP_SNP_uncor",  "FD_TC_bayP", "TP_TC_uncor","FD_WZA_bayP","TP_WZA_uncor", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP" ,"rep", "map")  )
Stabilising_FalsePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_FalsePositives_Uncor$variable <- factor(Stabilising_FalsePositives_Uncor$variable,
                                                    levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_uncor"),
                                                    labels = c( "SNP-based", "Top-Candidate", "WZA"))

Stabilising_FalsePositives_Baypass <- melt( Stabilising, id = c( "X" , "top" ,"TP_SNP_bayP", "FD_SNP_uncor",  "TP_TC_bayP", "FD_TC_uncor","TP_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor",  "rep", "map")  )
Stabilising_FalsePositives_Baypass$data <-  "BayPass"
Stabilising_FalsePositives_Baypass$variable <- factor(Stabilising_FalsePositives_Baypass$variable,
                                                      levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
                                                      labels = c( "SNP-based", "Top-Candidate", "WZA"))




Directional_FalsePositives_Uncor$selection <- "Directional Selection"
Stabilising_FalsePositives_Uncor$selection <- "Stabilising Selection"

FalsePositives_Uncor <- rbind(Directional_FalsePositives_Uncor, Stabilising_FalsePositives_Uncor)

Directional_FalsePositives_Baypass$selection <- "Directional Selection"
Stabilising_FalsePositives_Baypass$selection <- "Stabilising Selection"

TruePositives_Baypass <- rbind(Directional_FalsePositives_Baypass, Stabilising_FalsePositives_Baypass)



directional_falsePositives <- ggplot()+
  geom_line(data = Directional_FalsePositives_Uncor[Directional_FalsePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_FalsePositives_Baypass[Directional_FalsePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Directional Selection")+
  scale_colour_brewer('Test', palette = "Set2")+
  scale_fill_brewer('Test', palette = "Set2")+
  scale_y_continuous("False Discovery Rate (False Positives/Number of Top # Genes)", expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


stabilising_falsePositives <- ggplot()+
  geom_line(data = Stabilising_FalsePositives_Uncor[Stabilising_FalsePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Stabilising_FalsePositives_Baypass[Stabilising_FalsePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Stabilising Selection")+
  scale_colour_brewer('Test', palette = "Set2")+
  scale_fill_brewer('Test', palette = "Set2")+
  scale_y_continuous("False Discovery Rate (False Positives/Number of Top # Genes)", expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
    
  )

combinedPlot_FDR <- ggarrange(directional_falsePositives,stabilising_falsePositives, ncol= 2, common.legend = T, legend = "bottom", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_FalsePositives.pdf", height = 8, width = 7)
print(combinedPlot_FDR)
dev.off()







dir_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.csv")

