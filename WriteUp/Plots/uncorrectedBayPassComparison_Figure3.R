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
Directional_id_vars_uncor <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_uncor","TP_WZA_empR"))]

Directional_TruePositives_Uncor <- melt( Directional, id = Directional_id_vars_uncor  )
Directional_TruePositives_Uncor$data <- "Uncorrected p-values"
Directional_TruePositives_Uncor$variable <- factor(Directional_TruePositives_Uncor$variable,
                                                   levels = c( "TP_SNP_uncor",  "TP_WZA_empR"),
                                                   labels = c("B_TP_SNP_uncor","A_TP_WZA_empR"))

#,
#                                                   labels = c( "SNP-based",  "WZA"))

Directional_id_vars_BayPass <- names(Directional)[!(names(Directional)%in%c( "TP_SNP_bayP",  "TP_WZA_bayP"))]

Directional_TruePositives_Baypass <- melt( Directional, id = Directional_id_vars_BayPass )
Directional_TruePositives_Baypass$data <-  "BayPass"
Directional_TruePositives_Baypass$variable <- factor(Directional_TruePositives_Baypass$variable,
                                                     levels = c( "TP_SNP_bayP", "TP_WZA_bayP"),
                                                     labels = c("C_TP_SNP_bayP","D_TP_WZA_bayP"))

#,
#                                                     labels = c( "SNP-based", "Top-Candidate", "WZA"))



stab_bc_map <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.summary.csv")
stab_bc_map$map <- "BC Map"

stab_cline <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/cline_sampled.WZA.summary.csv")
stab_cline$map <- "Gradient"

stab_trunc <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/trunc_sampled.WZA.summary.csv")
stab_trunc$map <- "Truncated"

Stabilizing<- rbind(stab_bc_map, stab_cline, stab_trunc)

Stabilizing_id_vars_uncor <- names(Stabilizing)[!(names(Stabilizing)%in%c( "PCE_SNP_uncor", "PCE_WZA_empR"))]
Stabilizing_TruePositives_Uncor <- melt( Stabilizing, id = Stabilizing_id_vars_uncor  )
Stabilizing_TruePositives_Uncor$data <- "Uncorrected p-values"
Stabilizing_TruePositives_Uncor$variable <- factor(Stabilizing_TruePositives_Uncor$variable,
                                                   levels = c( "PCE_SNP_uncor",  "PCE_WZA_empR"),
                                                   labels = c("B_PCE_SNP_uncor","A_PCE_WZA_empR"))
#,
#                                                   labels = c( expression("Kendall's "*tau),   expression("WZA"*tau)))

Stabilizing_id_vars_BayPass <- names(Stabilizing)[!(names(Stabilizing)%in%c( "PCE_SNP_bayP", "PCE_WZA_bayP"))]
Stabilizing_TruePositives_Baypass <- melt( Stabilizing, id = Stabilizing_id_vars_BayPass)
Stabilizing_TruePositives_Baypass$data <-  "BayPass"
Stabilizing_TruePositives_Baypass$variable <- factor(Stabilizing_TruePositives_Baypass$variable,
                                                     levels = c( "PCE_SNP_bayP", "PCE_WZA_bayP"),
                                                     labels = c("D_PCE_SNP_bayP","C_PCE_WZA_bayP"))
#,
#                                                     labels = c( expression("BayPass"),  expression("WZA"[BP])))





Directional_TruePositives_Uncor$selection <- "Directional Selection"
Stabilizing_TruePositives_Uncor$selection <- "Stabilizing Selection"


Directional_TruePositives_Baypass$selection <- "Directional Selection"
Stabilizing_TruePositives_Baypass$selection <- "Stabilizing Selection"




paletteCB <- c("#D55E00",
               "#999999",
               "#CC79A7",
               "#009E73")
library(scales)

labs1 = c(expression("WZA"*tau),
          expression("Kendall's "*tau),
          expression("WZA"[BP]),
          expression("BayPass"))

directional_truePositives <- ggplot()+
  geom_line(data = Directional_TruePositives_Uncor[Directional_TruePositives_Uncor$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  geom_line(data = Directional_TruePositives_Baypass[Directional_TruePositives_Baypass$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  facet_grid(selection~ map)+
  theme_bw()+
  ggtitle("Directional Selection")+
  scale_color_manual("Method",
                     values = paletteCB,
                     labels = labs1)+
  
  scale_y_continuous(expression("Proportion of true positives detected"), expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    strip.text.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.text.align = 0,   
    plot.title = element_blank()
    #plot.title = element_text(hjust = 0.5, size = 15)
  )

stablising_truePositives <- ggplot()+
  geom_line(data = Stabilizing_TruePositives_Uncor[Stabilizing_TruePositives_Uncor$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  geom_line(data = Stabilizing_TruePositives_Baypass[Stabilizing_TruePositives_Baypass$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  facet_grid(selection~ map)+
  theme_bw()+
  ggtitle("Stablising Selection")+
  scale_color_manual("Method",
                      values = paletteCB,
                     labels = labs1)+
                     
  scale_y_continuous(expression("Proportion of Cov["*italic("phen, env")*"] explained"), expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    strip.text.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.text.align = 0,
    plot.title = element_blank()
#    plot.title = element_text(hjust = 0.5, size = 15)
  )

combinedPlot <- ggarrange(annotate_figure(directional_truePositives, top = "Directional Selection"), 
                          annotate_figure(stablising_truePositives, top = "Stabilizing Selection"), 
                          ncol= 1, nrow =2, common.legend = T,
                          legend = "right", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_TruePositives.pdf", height = 6, width = 6)
print(combinedPlot)
dev.off()






















### Now do the false positives


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




Stabilizing_FalsePositives_Uncor_ID_Vars <- names(Stabilizing)[!(names(Stabilizing)%in%c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"))]
Stabilizing_FalsePositives_Uncor <- melt( Stabilizing, id = Stabilizing_FalsePositives_Uncor_ID_Vars  )
Stabilizing_FalsePositives_Uncor$data <- "Uncorrected p-values"
Stabilizing_FalsePositives_Uncor$variable <- factor(Stabilizing_FalsePositives_Uncor$variable,
                                                    levels = c( "FD_SNP_uncor", "FD_TC_uncor", "FD_WZA_empR"),
                                                    labels = c( "SNP-based", "Top-Candidate", "WZA"))

Stabilizing_FalsePositives_Uncor_ID_Vars <- names(Stabilizing)[!(names(Stabilizing)%in%c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"))]
Stabilizing_FalsePositives_Baypass <- melt( Stabilizing, id = c( "X" , "top" ,"TP_SNP_bayP", "FD_SNP_uncor",  "TP_TC_bayP", "FD_TC_uncor","TP_WZA_bayP","FD_WZA_uncor",   "TP_SNP_uncor",   "TP_TC_uncor",  "TP_WZA_uncor", "TP_WZA_empR",  "rep", "map")  )
Stabilizing_FalsePositives_Baypass$data <-  "BayPass"
Stabilizing_FalsePositives_Baypass$variable <- factor(Stabilizing_FalsePositives_Baypass$variable,
                                                      levels = c( "FD_SNP_bayP", "FD_TC_bayP", "FD_WZA_bayP"),
                                                      labels = c( "SNP-based", "Top-Candidate", "WZA"))




Directional_FalsePositives_Uncor$selection <- "Directional Selection"
Stabilizing_FalsePositives_Uncor$selection <- "Stabilizing Selection"


Directional_FalsePositives_Baypass$selection <- "Directional Selection"
Stabilizing_FalsePositives_Baypass$selection <- "Stabilizing Selection"




directional_falsePositives <- ggplot()+
  geom_line(data = Directional_FalsePositives_Uncor[Directional_FalsePositives_Uncor$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 0.9)+
  geom_line(data = Directional_FalsePositives_Baypass[Directional_FalsePositives_Baypass$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 0.9)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Directional Selection")+
  scale_color_manual("Method", values = paletteCB)+
  scale_y_continuous("False Discovery Rate\n(False Positives/Number of Number of Genes in Top Set)", expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 30),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


stabilising_falsePositives <- ggplot()+
  geom_line(data = Stabilizing_FalsePositives_Uncor[Stabilizing_FalsePositives_Uncor$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  geom_line(data = Stabilizing_FalsePositives_Baypass[Stabilizing_FalsePositives_Baypass$rep == "mean",], aes( x = top/1000, y = value, col = variable), lwd = 1)+
  facet_grid(data ~ map)+
  theme_bw()+
  ggtitle("Stabilizing Selection")+
  scale_color_manual("Method", values = paletteCB)+
  scale_y_continuous("False discovery rate \n(False Positives/Number of Number of Genes in Top Set)", expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 30),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
    
  )

combinedPlot_FDR <- ggarrange(directional_falsePositives,stabilising_falsePositives, ncol= 2, common.legend = T, legend = "right", labels = "AUTO")

pdf("~/work/GEA/simulations/Plots/UncorrectedBayPassComparison_FalsePositives.pdf", height = 4, width = 10)
print(combinedPlot_FDR)
dev.off()



## Make a plot of ProbGen poster
probGenData_Uncor <- Stabilizing_TruePositives_Uncor[Stabilizing_TruePositives_Uncor$rep == "mean",]
probGenData_Uncor <- probGenData_Uncor[probGenData_Uncor$map == "BC Map",]
probGenData_Uncor <- probGenData_Uncor[probGenData_Uncor$variable != "Top-Candidate",]

probGenData_Baypass <- Stabilizing_TruePositives_Baypass[Stabilizing_TruePositives_Baypass$rep == "mean",]
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
  geom_line(data = probGenData_Uncor, aes( x = top/1000, y = value, col = variable), lwd = 1)+
  geom_line(data = probGenData_Baypass, aes( x = top/1000, y = value, col = variable), lwd = 1)+
  facet_grid( ~ data, drop = T)+
  theme_bw()+
  ggtitle("Stablizing Selection")+
  scale_color_manual("Method", values = c(paletteCB[2],paletteCB[4]))+
  scale_y_continuous(expression("Proportion of local adaptation explained"), expand = c(0,0))+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )





dir_BC_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.csv")
dir_BC_wza$map <- "BC map"
dir_cline_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/cline_sampled.WZA.csv")
dir_cline_wza$map <- "Gradient map"
dir_trunc_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.csv")
dir_trunc_wza$map <- "Truncated map"
dir <- rbind(dir_BC_wza,  dir_cline_wza, dir_trunc_wza )
dir$selection <- "Directional Selection"

stab_BC_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.csv")
stab_BC_wza$map <- "BC map"
stab_cline_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/cline_sampled.WZA.csv")
stab_cline_wza$map <- "Gradient map"
stab_trunc_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/trunc_sampled.WZA.csv")
stab_trunc_wza$map <- "Truncated map"
stab <- rbind(stab_BC_wza,  stab_cline_wza, stab_trunc_wza )
stab$selection <- "Stabilizing Selection"

sel <- rbind( dir, stab )

Z_score_distribution <- ggplot(data = sel, aes( x = empR_Z, fill = LA > 0.005))+
  geom_density(alpha = 0.75)+
  scale_y_continuous("Density", limits = c(0,0.5))+
  scale_x_continuous(expression(italic(Z[W])))+
  facet_grid(selection~map)+
  scale_fill_manual("Locally adaptive gene", values = c("black", "lightblue"))+
  theme_bw()+
  theme(
    strip.text = element_text(face= "italic", size = 10),
    strip.background = element_blank()
  )

pdf("~/work/GEA/simulations/Plots/Z_score_distribution_selection.pdf", width = 10, height = 7)
print(Z_score_distribution)
dev.off()
