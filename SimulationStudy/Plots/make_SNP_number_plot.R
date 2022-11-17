rm(list = ls())

library(ggplot2)

bc_strong_d10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20_SNP_number.csv")
bc_strong_d10_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20")
bc_strong_d10_comb <- rbind(bc_strong_d10, bc_strong_d10_all[bc_strong_d10_all$stat=="Z_geno_k_tau_p_value",])
bc_strong_d10_comb$demes <- "10 Demes"

bc_strong_d20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20_SNP_number.csv")
bc_strong_d20_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20")
bc_strong_d20_comb <- rbind(bc_strong_d20, bc_strong_d20_all[bc_strong_d20_all$stat=="Z_geno_k_tau_p_value",])
bc_strong_d20_comb$demes <- "20 Demes"

bc_strong_d40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20_SNP_number.csv")
bc_strong_d40_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20")
bc_strong_d40_comb <- rbind(bc_strong_d40, bc_strong_d40_all[bc_strong_d40_all$stat=="Z_geno_k_tau_p_value",])
bc_strong_d40_comb$demes <- "40 Demes"

bc_strong <- rbind( bc_strong_d10_comb, 
                    bc_strong_d20_comb,
                    bc_strong_d40_comb)
bc_strong$map <- "BC Map"




cline_strong_d10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20_SNP_number.csv")
cline_strong_d10_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20")
cline_strong_d10_comb <- rbind(cline_strong_d10, cline_strong_d10_all[cline_strong_d10_all$stat=="Z_geno_k_tau_p_value",])
cline_strong_d10_comb$demes <- "10 Demes"

cline_strong_d20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20_SNP_number.csv")
cline_strong_d20_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20")
cline_strong_d20_comb <- rbind(cline_strong_d20, cline_strong_d20_all[cline_strong_d20_all$stat=="Z_geno_k_tau_p_value",])
cline_strong_d20_comb$demes <- "20 Demes"

cline_strong_d40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20_SNP_number.csv")
cline_strong_d40_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20")
cline_strong_d40_comb <- rbind(cline_strong_d40, cline_strong_d40_all[cline_strong_d40_all$stat=="Z_geno_k_tau_p_value",])
cline_strong_d40_comb$demes <- "40 Demes"

cline_strong <- rbind( cline_strong_d10_comb, 
                    cline_strong_d20_comb,
                    cline_strong_d40_comb)
cline_strong$map <- "Gradient Map"



trunc_strong_d10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20_SNP_number.csv")
trunc_strong_d10_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20")
trunc_strong_d10_comb <- rbind(trunc_strong_d10, trunc_strong_d10_all[trunc_strong_d10_all$stat=="Z_geno_k_tau_p_value",])
trunc_strong_d10_comb$demes <- "10 Demes"

trunc_strong_d20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20_SNP_number.csv")
trunc_strong_d20_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20")
trunc_strong_d20_comb <- rbind(trunc_strong_d20, trunc_strong_d20_all[trunc_strong_d20_all$stat=="Z_geno_k_tau_p_value",])
trunc_strong_d20_comb$demes <- "20 Demes"

trunc_strong_d40 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20_SNP_number.csv")
trunc_strong_d40_all <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20")
trunc_strong_d40_comb <- rbind(trunc_strong_d40, trunc_strong_d40_all[trunc_strong_d40_all$stat=="Z_geno_k_tau_p_value",])
trunc_strong_d40_comb$demes <- "40 Demes"

trunc_strong <- rbind( trunc_strong_d10_comb, 
                       trunc_strong_d20_comb,
                       trunc_strong_d40_comb)
trunc_strong$map <- "Truncated Map"

all <- rbind(bc_strong,cline_strong,trunc_strong)


parse.labels <- function(x) parse(text = x)

colorBlindGrey8 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

all <- all[all$stat %in% c("geno_k_tau_p_value",
                           "Z_10snp",
                           "Z_2snp",
                           "Z_4snp",
                           "Z_8snp",
                           "Z_geno_k_tau_p_value"),]

all$statistic <- factor(all$stat,
                        levels = c("geno_k_tau_p_value",
                                   "Z_2snp",
                                   "Z_4snp",
                                   "Z_8snp",
                                   "Z_10snp",
                                   "Z_geno_k_tau_p_value"),
                        labels = c(expression("Kendall's "*tau*" (Single SNP)"),
                                   expression("WZA"[tau]*" (2 SNPs)"),
                                   expression("WZA"[tau]*" (4 SNPs)"),
                                   expression("WZA"[tau]*" (8 SNPs)"),
                                   expression("WZA"[tau]*" (10 SNPs)"),
                                   expression("WZA"[tau]*" (10,000bp)")))


snp_number_plot <- ggplot(data = all, 
                                                   aes(x = statistic,
                                                       y = auc_pr))+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("~/work/GEA/simulations/SimulationStudy/Plots/SNP_number_plot.pdf",
       plot = snp_number_plot,
       width = 9.08,
       height = 5.71)
