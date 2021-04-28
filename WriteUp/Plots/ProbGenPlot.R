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
Directional_id_vars_uncor <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_uncor", "TP_WZA_empR"))]

Directional_TruePositives_Uncor <- melt( Directional, id = Directional_id_vars_uncor  )
Directional_TruePositives_Uncor$data <- "Uncorrected p-values"
Directional_TruePositives_Uncor$variable <- factor(Directional_TruePositives_Uncor$variable,
                                                   levels = c( "TP_SNP_uncor", "TP_WZA_empR"),
                                                   labels = c( "SNP-based", "WZA"))

Directional_id_vars_BayPass <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_bayP", "TP_WZA_bayP"))]

Directional_TruePositives_Baypass <- melt( Directional, id = Directional_id_vars_BayPass )
Directional_TruePositives_Baypass$data <-  "BayPass"
Directional_TruePositives_Baypass$variable <- factor(Directional_TruePositives_Baypass$variable,
                                                     levels = c( "TP_SNP_bayP", "TP_WZA_bayP"),
                                                     labels = c( "SNP-based",  "WZA"))



paletteCB <- c("#1d457f",
               "#cc5c76",
               "#f9ad2a",
               "#0571b0")

probGenPlot_1 <- ggplot()+
  geom_line(data = Directional_TruePositives_Uncor[Directional_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_TruePositives_Baypass[Directional_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid( data ~ map )+
  theme_bw()+
  ggtitle("Correlation between environment and selection = 1.0")+
  scale_color_manual("Test", values = paletteCB)+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 10),
    #strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 13)
  )
png("~/work/GEA/simulations/Plots/ProbGenPoster_Plot1.png", width = 8, height = 5, units = "in", res = 300)
print( probGenPlot_1 )
dev.off()


dir_bc_map_cor_0.5 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled_cor0.5.WZA.summary.csv")
dir_bc_map_cor_0.5 <- dir_bc_map_cor_0.5[dir_bc_map_cor_0.5$rep == "mean",]
dir_bc_map_cor_0.5$map <- "BC Map"


Stabilising<- rbind(dir_bc_map_cor_0.5)

Stabilising_id_vars_uncor <- names(Stabilising)[!(names(Stabilising)%in%c( "PCE_SNP_uncor",  "PCE_WZA_empR"))]
Stabilising_TruePositives_Uncor <- melt( Stabilising, id = Stabilising_id_vars_uncor  )
Stabilising_TruePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_TruePositives_Uncor$variable <- factor(Stabilising_TruePositives_Uncor$variable,
                                                   levels = c( "PCE_SNP_uncor", "PCE_WZA_empR"),
                                                   labels = c( "SNP-based",  "WZA"))

Stabilising_id_vars_BayPass <- names(Stabilising)[!(names(Stabilising)%in%c( "PCE_SNP_bayP",  "PCE_WZA_bayP"))]
Stabilising_TruePositives_Baypass <- melt( Stabilising, id = Stabilising_id_vars_BayPass  )
Stabilising_TruePositives_Baypass$data <-  "BayPass"
Stabilising_TruePositives_Baypass$variable <- factor(Stabilising_TruePositives_Baypass$variable,
                                                     levels = c( "PCE_SNP_bayP", "PCE_WZA_bayP"),
                                                     labels = c( "SNP-based", "WZA"))


envCor_Plot <- ggplot()+
  geom_line(data = Stabilising_TruePositives_Uncor, aes( x = top, y = value, col = variable), lwd = 1)+
  geom_line(data = Stabilising_TruePositives_Baypass, aes( x = top, y = value, col = variable), lwd = 1)+
  facet_grid(data ~ map, drop = T)+
  theme_bw()+
  ggtitle("Correlation between\nenvironment and selection = 0.5")+
  scale_color_manual("Test", values = c(paletteCB[4],paletteCB[2]))+
  scale_y_continuous(expression("Proportion of local adaptation explained"), expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size =13)
  )


library(ggpubr)

png("~/work/GEA/simulations/Plots/ProbGenPoster_Plot2.png", width = 9, height = 5, units = "in", res = 300)
ggarrange(probGenPlot_1, envCor_Plot, ncol = 2, widths = c(2.2,1),common.legend = F, align = "h", legend = "bottom", labels = "AUTO" )
dev.off()
