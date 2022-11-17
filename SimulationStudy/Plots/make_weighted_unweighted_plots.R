rm(list = ls())

bc_map_strong_10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20")
bc_map_strong_10$n <-"10 Demes"
bc_map_strong_20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20")
bc_map_strong_20$n <- "20 Demes"
bc_map_strong_40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20")
bc_map_strong_40$n <- "40 Demes"
bc_map_strong <- rbind(bc_map_strong_10,bc_map_strong_20,bc_map_strong_40)
bc_map_strong$map<-"BC Map"

cline_strong_10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20")
cline_strong_10$n <-"10 Demes"
cline_strong_20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20")
cline_strong_20$n <-"20 Demes"
cline_strong_40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20")
cline_strong_40$n <-"40 Demes"
cline_strong <- rbind(cline_strong_10,cline_strong_20,cline_strong_40)
cline_strong$map<-"Gradient Map"

trunc_strong_10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20")
trunc_strong_10$n <-"10 Demes"
trunc_strong_20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20")
trunc_strong_20$n <- "20 Demes"
trunc_strong_40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20")
trunc_strong_40$n <- "40 Demes"
trunc_strong <- rbind(trunc_strong_10,trunc_strong_20,trunc_strong_40)
trunc_strong$map<-"Truncated Map"

vanilla_WZA <- rbind( bc_map_strong, cline_strong, trunc_strong)
vanilla_WZA <- vanilla_WZA[vanilla_WZA$stat == "Z_geno_k_tau_p_value", ]

## Supplementary Stats

bc_map_strong_10_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20_supp.csv")
bc_map_strong_10_s$n <-"10 Demes"
bc_map_strong_20_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20_supp.csv")
bc_map_strong_20_s$n <- "20 Demes"
bc_map_strong_40_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20_supp.csv")
bc_map_strong_40_s$n <- "40 Demes"
bc_map_strong_s <- rbind(bc_map_strong_10_s,bc_map_strong_20_s,bc_map_strong_40_s)
bc_map_strong_s$map<-"BC Map"

cline_strong_10_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20_supp.csv")
cline_strong_10_s$n <-"10 Demes"
cline_strong_20_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20_supp.csv")
cline_strong_20_s$n <-"20 Demes"
cline_strong_40_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20_supp.csv")
cline_strong_40_s$n <-"40 Demes"
cline_strong_s <- rbind(cline_strong_10_s,cline_strong_20_s,cline_strong_40_s)
cline_strong_s$map<-"Gradient Map"

trunc_strong_10_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20_supp.csv")
trunc_strong_10_s$n <-"10 Demes"
trunc_strong_20_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20_supp.csv")
trunc_strong_20_s$n <- "20 Demes"
trunc_strong_40_s <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20_supp.csv")
trunc_strong_40_s$n <- "40 Demes"
trunc_strong_s <- rbind(trunc_strong_10_s,trunc_strong_20_s,trunc_strong_40_s)
trunc_strong_s$map<-"Truncated Map"

supp_WZA <- rbind( bc_map_strong_s, cline_strong_s, trunc_strong_s)


# Script to make the figure showing the weighted/unweighted WZA

supp_WZA_unweighted <- supp_WZA[supp_WZA$stat == "Z_geno_k_tau_p_value_unweighted",]
WZA_unweighted_dat <- rbind(vanilla_WZA, supp_WZA_unweighted)

library(ggplot2)

labels_main_figs <- c(expression("WZA"[tau]*" - Weighted"),
                      expression("WZA"[tau]*" - Unweighted"))

main_stats <- c("Z_geno_k_tau_p_value", "Z_geno_k_tau_p_value_unweighted")

WZA_unweighted_dat$statistic  <- factor(WZA_unweighted_dat$stat,
                         levels = main_stats,
                         labels = labels_main_figs)



parse.labels <- function(x) parse(text = x)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "white")

weighted_unweighted  <- ggplot(data = WZA_unweighted_dat, 
                                                   aes(x = statistic,
                                                       y = auc_pr))+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~n)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =c(colorBlindGrey8[6],"black"))+
  scale_fill_manual(labels = parse.labels,
                    values =c(colorBlindGrey8[6],colorBlindGrey8[9]))+
  theme_bw()+

  theme(axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))


ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Weighted_Unweighted_Comparison.pdf",
       plot = weighted_unweighted,
       width = 9.08,
       height = 5.71)

# Script to show the benefit to empiricalal p-values to parametric p-values


supp_WZA_parametric <- supp_WZA[supp_WZA$stat == "Z_parametri",]
WZA_parametric_dat <- rbind(vanilla_WZA, supp_WZA_parametric)

labels_main_figs <- c(expression("WZA"[tau]*" - Empirical "*italic(p)*"-values"),
                      expression("WZA"[tau]*" - Parametric "*italic(p)*"-values"))

main_stats <- c("Z_geno_k_tau_p_value", "Z_parametri")

WZA_parametric_dat$statistic  <- factor(WZA_parametric_dat$stat,
                                        levels = main_stats,
                                        labels = labels_main_figs)


parse.labels <- function(x) parse(text = x)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "white")

empirical_parametric  <- ggplot(data = WZA_parametric_dat, 
                               aes(x = statistic,
                                   y = auc_pr))+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~n, scales = "free_y")+
#  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =c(colorBlindGrey8[6],"black"))+
  scale_fill_manual(labels = parse.labels,
                    values =c(colorBlindGrey8[6],colorBlindGrey8[9]))+
  theme_bw()+
  
  theme(axis.text.x = element_text(color = "black",
                                   angle = 80, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Empirical_Parametric_Comparison.pdf",
       plot = empirical_parametric,
       width = 9.08,
       height = 7)
