rm(list = ls())

## Read in the necessary packages
library(reshape2)
library(ggplot2)
library(ggpubr)

dir_bc_map <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.summary.csv")
dir_bc_map$map <- "BC Map"
<<<<<<< HEAD

=======
  
>>>>>>> origin/master
dir_cline <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/cline_sampled.WZA.summary.csv")
dir_cline$map <- "Gradient"

dir_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.summary.csv")
dir_trunc$map <- "Truncated"

Directional<- rbind(dir_bc_map, dir_cline, dir_trunc)
<<<<<<< HEAD
Directional_id_vars_uncor <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_empR"))]

Directional_TruePositives_Uncor <- melt( Directional, id = Directional_id_vars_uncor  )
Directional_TruePositives_Uncor$data <- "Uncorrected p-values"
Directional_TruePositives_Uncor$variable <- factor(Directional_TruePositives_Uncor$variable,
                                                   levels = c( "TP_SNP_uncor", "TP_TC_uncor", "TP_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA"))

Directional_id_vars_BayPass <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_bayP", "TP_TC_bayP", "TP_WZA_bayP"))]

Directional_TruePositives_Baypass <- melt( Directional, id = Directional_id_vars_BayPass )
Directional_TruePositives_Baypass$data <-  "BayPass"
Directional_TruePositives_Baypass$variable <- factor(Directional_TruePositives_Baypass$variable,
                                                     levels = c( "TP_SNP_bayP", "TP_TC_bayP", "TP_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA"))
=======

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
>>>>>>> origin/master



stab_bc_map <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.summary.csv")
stab_bc_map$map <- "BC Map"

stab_cline <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/cline_sampled.WZA.summary.csv")
stab_cline$map <- "Gradient"

stab_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/trunc_sampled.WZA.summary.csv")
stab_trunc$map <- "Truncated"

Stabilising<- rbind(stab_bc_map, stab_cline, stab_trunc)

<<<<<<< HEAD
Stabilising_id_vars_uncor <- names(Stabilising)[!(names(Stabilising)%in%c( "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR"))]
Stabilising_TruePositives_Uncor <- melt( Stabilising, id = Stabilising_id_vars_uncor  )
Stabilising_TruePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_TruePositives_Uncor$variable <- factor(Stabilising_TruePositives_Uncor$variable,
                                                   levels = c( "PCE_SNP_uncor", "PCE_TC_uncor", "PCE_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate",  "WZA"))

Stabilising_id_vars_BayPass <- names(Stabilising)[!(names(Stabilising)%in%c( "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP"))]
Stabilising_TruePositives_Baypass <- melt( Stabilising, id = c( "X" , "top" ,"FD_SNP_bayP", "FD_SNP_uncor",  "FD_TC_bayP", "FD_TC_uncor","FD_WZA_bayP","FD_WZA_uncor","FD_WZA_empR",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR", "TP_SNP_bayP", "TP_TC_bayP", "TP_WZA_bayP","rep", "map")  )
Stabilising_TruePositives_Baypass$data <-  "BayPass"
Stabilising_TruePositives_Baypass$variable <- factor(Stabilising_TruePositives_Baypass$variable,
                                                     levels = c( "PCE_SNP_bayP", "PCE_TC_bayP", "PCE_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA"))
=======
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
>>>>>>> origin/master





Directional_TruePositives_Uncor$selection <- "Directional Selection"
Stabilising_TruePositives_Uncor$selection <- "Stabilising Selection"

<<<<<<< HEAD
=======
TruePositives_Uncor <- rbind(Directional_TruePositives_Uncor, Stabilising_TruePositives_Uncor)
>>>>>>> origin/master

Directional_TruePositives_Baypass$selection <- "Directional Selection"
Stabilising_TruePositives_Baypass$selection <- "Stabilising Selection"

<<<<<<< HEAD


paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")

=======
TruePositives_Baypass <- rbind(Directional_TruePositives_Baypass, Stabilising_TruePositives_Baypass)

paletteCB <- c("#ca0020",
               "#f4a582",
               "#92c5de",
               "#0571b0")
>>>>>>> origin/master

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
<<<<<<< HEAD
  geom_line(data = Stabilising_TruePositives_Uncor[Stabilising_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 1)+
  geom_line(data = Stabilising_TruePositives_Baypass[Stabilising_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 1)+
=======
  geom_line(data = Stabilising_TruePositives_Uncor[Stabilising_TruePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Stabilising_TruePositives_Baypass[Stabilising_TruePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
>>>>>>> origin/master
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Stablising Selection")+
  scale_color_manual("Test", values = paletteCB)+
<<<<<<< HEAD
  scale_y_continuous(expression("Proportion of Cov["*italic("phen, env")*"] Explained"), expand = c(0,0))+
=======
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
>>>>>>> origin/master
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
<<<<<<< HEAD
    plot.title = element_text(hjust = 0.5, size = 15)
=======
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
    
>>>>>>> origin/master
  )

combinedPlot <- ggarrange(directional_truePositives,stablising_truePositives, ncol= 2, common.legend = T, legend = "right", labels = "AUTO")

<<<<<<< HEAD
pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_TruePositives.pdf", height = 4, width = 10)
=======
pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_TruePositives.pdf", height = 9.5, width = 9)
>>>>>>> origin/master
print(combinedPlot)
dev.off()



### Now do the false positives


<<<<<<< HEAD
Directional_FalsePositives_Uncor_ID_Vars <- names(Directional)[!(names(Directional)%in%c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"))]
Directional_FalsePositives_Uncor <- melt( Directional, id = Directional_FalsePositives_Uncor_ID_Vars  )
Directional_FalsePositives_Uncor$data <- "Uncorrected p-values"
Directional_FalsePositives_Uncor$variable <- factor(Directional_FalsePositives_Uncor$variable,
                                                   levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA"))

Directional_FalsePositives_BayPass_ID_Vars <- names(Directional)[!(names(Directional)%in%c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"))]
Directional_FalsePositives_Baypass <- melt( Directional, id = Directional_FalsePositives_BayPass_ID_Vars  )
Directional_FalsePositives_Baypass$data <-  "BayPass"
Directional_FalsePositives_Baypass$variable <- factor(Directional_FalsePositives_Baypass$variable,
                                                     levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA"))
=======

Directional_FalsePositives_Uncor <- melt( Directional, id = c( "X" , "top" ,"FD_SNP_bayP", "TP_SNP_uncor",  "FD_TC_bayP", "TP_TC_uncor","FD_WZA_bayP","TP_WZA_uncor", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP", "TP_WZA_empR" ,"rep", "map")  )
Directional_FalsePositives_Uncor$data <- "Uncorrected p-values"
Directional_FalsePositives_Uncor$variable <- factor(Directional_FalsePositives_Uncor$variable,
                                                   levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_uncor", "FD_WZA_empR"),
                                                   labels = c( "SNP-based", "Top-Candidate", "WZA - p-values", "WZA - empRho"))

Directional_FalsePositives_Baypass <- melt( Directional, id = c( "X" , "top" ,"TP_SNP_bayP", "FD_SNP_uncor",  "TP_TC_bayP", "FD_TC_uncor","TP_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR", "rep", "map")  )
Directional_FalsePositives_Baypass$data <-  "BayPass"
Directional_FalsePositives_Baypass$variable <- factor(Directional_FalsePositives_Baypass$variable,
                                                     levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
                                                     labels = c( "SNP-based", "Top-Candidate", "WZA - p-values"))
>>>>>>> origin/master




<<<<<<< HEAD
Stabilising_FalsePositives_Uncor_ID_Vars <- names(Stabilising)[!(names(Stabilising)%in%c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"))]
Stabilising_FalsePositives_Uncor <- melt( Stabilising, id = Stabilising_FalsePositives_Uncor_ID_Vars  )
Stabilising_FalsePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_FalsePositives_Uncor$variable <- factor(Stabilising_FalsePositives_Uncor$variable,
                                                    levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"),
                                                    labels = c( "SNP-based", "Top-Candidate", "WZA"))

Stabilising_FalsePositives_Uncor_ID_Vars <- names(Stabilising)[!(names(Stabilising)%in%c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"))]
=======
Stabilising_FalsePositives_Uncor <- melt( Stabilising, id = c( "X" , "top" ,"FD_SNP_bayP", "TP_SNP_uncor",  "FD_TC_bayP", "TP_TC_uncor","FD_WZA_bayP","TP_WZA_uncor", "TP_SNP_bayP", "TP_TC_bayP","TP_WZA_bayP", "TP_WZA_empR","rep", "map")  )
Stabilising_FalsePositives_Uncor$data <- "Uncorrected p-values"
Stabilising_FalsePositives_Uncor$variable <- factor(Stabilising_FalsePositives_Uncor$variable,
                                                    levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_uncor", "FD_WZA_empR"),
                                                    labels = c( "SNP-based", "Top-Candidate", "WZA - p-values", "WZA - empRho"))

>>>>>>> origin/master
Stabilising_FalsePositives_Baypass <- melt( Stabilising, id = c( "X" , "top" ,"TP_SNP_bayP", "FD_SNP_uncor",  "TP_TC_bayP", "FD_TC_uncor","TP_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR",  "rep", "map")  )
Stabilising_FalsePositives_Baypass$data <-  "BayPass"
Stabilising_FalsePositives_Baypass$variable <- factor(Stabilising_FalsePositives_Baypass$variable,
                                                      levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
<<<<<<< HEAD
                                                      labels = c( "SNP-based", "Top-Candidate", "WZA"))
=======
                                                      labels = c( "SNP-based", "Top-Candidate", "WZA - p-values"))
>>>>>>> origin/master




Directional_FalsePositives_Uncor$selection <- "Directional Selection"
Stabilising_FalsePositives_Uncor$selection <- "Stabilising Selection"

<<<<<<< HEAD
=======
FalsePositives_Uncor <- rbind(Directional_FalsePositives_Uncor, Stabilising_FalsePositives_Uncor)
>>>>>>> origin/master

Directional_FalsePositives_Baypass$selection <- "Directional Selection"
Stabilising_FalsePositives_Baypass$selection <- "Stabilising Selection"

<<<<<<< HEAD
=======
FalsePositives_Baypass <- rbind(Directional_FalsePositives_Baypass, Stabilising_FalsePositives_Baypass)
>>>>>>> origin/master



directional_falsePositives <- ggplot()+
  geom_line(data = Directional_FalsePositives_Uncor[Directional_FalsePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_FalsePositives_Baypass[Directional_FalsePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Directional Selection")+
  scale_color_manual("Test", values = paletteCB)+
<<<<<<< HEAD
  scale_y_continuous("False Discovery Rate\n(False Positives/Number of Top # Genes)", expand = c(0,0))+
=======
  scale_y_continuous("False Discovery Rate (False Positives/Number of Top # Genes)", expand = c(0,0))+
>>>>>>> origin/master
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
<<<<<<< HEAD
  geom_line(data = Stabilising_FalsePositives_Uncor[Stabilising_FalsePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 1)+
  geom_line(data = Stabilising_FalsePositives_Baypass[Stabilising_FalsePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 1)+
=======
  geom_line(data = Stabilising_FalsePositives_Uncor[Stabilising_FalsePositives_Uncor$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Stabilising_FalsePositives_Baypass[Stabilising_FalsePositives_Baypass$rep == "mean",], aes( x = top, y = value, col = variable), lwd = 0.9)+
>>>>>>> origin/master
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Stabilising Selection")+
  scale_color_manual("Test", values = paletteCB)+
<<<<<<< HEAD
  scale_y_continuous("False Discovery Rate \n(False Positives/Number of Top # Genes)", expand = c(0,0))+
=======
  scale_y_continuous("False Discovery Rate (False Positives/Number of Top # Genes)", expand = c(0,0))+
>>>>>>> origin/master
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

<<<<<<< HEAD
combinedPlot_FDR <- ggarrange(directional_falsePositives,stabilising_falsePositives, ncol= 2, common.legend = T, legend = "right", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_FalsePositives.pdf", height = 4, width = 10)
=======
combinedPlot_FDR <- ggarrange(directional_falsePositives,stabilising_falsePositives, ncol= 2, common.legend = T, legend = "bottom", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_FalsePositives.pdf", height = 9.5, width = 9)
>>>>>>> origin/master
print(combinedPlot_FDR)
dev.off()



<<<<<<< HEAD
## Make a plot of ProbGen poster
probGenData_Uncor <- Stabilising_TruePositives_Uncor[Stabilising_TruePositives_Uncor$rep == "mean",]
probGenData_Uncor <- probGenData_Uncor[probGenData_Uncor$map == "BC Map",]
probGenData_Uncor <- probGenData_Uncor[probGenData_Uncor$variable != "Top-Candidate",]

probGenData_Baypass <- Stabilising_TruePositives_Baypass[Stabilising_TruePositives_Baypass$rep == "mean",]
probGenData_Baypass <- probGenData_Baypass[probGenData_Baypass$map == "BC Map",]
probGenData_Baypass <- probGenData_Baypass[probGenData_Baypass$variable != "Top-Candidate",]


probGenData_Uncor$variable <- factor(probGenData_Uncor$variable,
                                       levels = c( "SNP-based", "WZA"),
                                       labels = c( "SNP-based", "WZA"))

probGenData_Uncor$data <- factor(probGenData_Uncor$data,
                                   levels = c( "Uncorrected p-values"),
                                   labels = c( "Uncorrected p-values"))


probGenData_Baypass$variable <- factor(probGenData_Baypass$variable,
                                                     levels = c( "SNP-based", "WZA"),
                                                     labels = c( "SNP-based", "WZA"))

probGenData_Baypass$data <- factor(probGenData_Baypass$data,
                                       levels = c( "BayPass"),
                                       labels = c( "BayPass"))

stablising_truePositives <- ggplot()+
  geom_line(data = probGenData_Uncor, aes( x = top, y = value, col = variable), lwd = 1)+
  geom_line(data = probGenData_Baypass, aes( x = top, y = value, col = variable), lwd = 1)+
  facet_grid( ~ data, drop = T)+
  theme_bw()+
  ggtitle("Stablising Selection")+
  scale_color_manual("Test", values = c(paletteCB[2],paletteCB[4]))+
  scale_y_continuous(expression("Proportion of local adaptation explained"), expand = c(0,0))+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )
=======

>>>>>>> origin/master
