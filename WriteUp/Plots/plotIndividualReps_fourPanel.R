rm(list = ls())

library(ggplot2)

BC_snp_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/BC_Map/sample40/1_0.003_1.12Loci.directionalSelection.sample.csv.csv")
BC_WZA_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/BC_Map/WZA_BC_Map_sample40.csv")
BC_summary_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/BC_Map/WZA_BC_Map_sample40.summary.csv")


BC_allSNPs <- ggplot(data = BC_snp_file[BC_snp_file$maf > 0.05,], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  geom_point()+
  geom_point(data = BC_snp_file[(BC_snp_file$maf > 0.05)&(BC_snp_file$LA > 0.005),], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(italic(-log[10]*"(Kendall's "*tau*" p-value)")))+
  theme_bw()

BC_WZA <- ggplot(data = BC_WZA_file, aes(x = position/1e6, y = empR_Z, col = LA> 0.005))+
  geom_point()+
  geom_point(data = BC_WZA_file[(BC_WZA_file$LA > 0.005),], aes(x = position/1e6, y = empR_Z, col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(Z[W]))+  
  theme_bw()


BC_TC <- ggplot(data = BC_WZA_file, aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005))+
  geom_point()+
  geom_point(data = BC_WZA_file[(BC_WZA_file$LA > 0.005),], aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous("Top-candidate\nIndex")+
  theme_bw()


paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")


BC_powerCurve <- ggplot(data = BC_summary_file[BC_summary_file$rep == "mean",], aes(x = top, y = TP_WZA_empR))+
  geom_line(aes(col = "WZA"),lwd = 2)+
  geom_line(data = BC_summary_file[BC_summary_file$rep == "mean",], aes(x = top, y = TP_TC_uncor, col = "Top-Candidate"),lwd = 2)+
  geom_line(data = BC_summary_file[BC_summary_file$rep == "mean",], aes(x = top, y = TP_SNP_uncor, col = "SNP-based"),lwd = 2)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  scale_x_continuous("Number of Genes in Top Set")+
  scale_color_manual("Test", values = paletteCB)+
  theme_bw()

library(ggpubr)




BC_panel_plot <- ggarrange( ggarrange(BC_allSNPs, labels = c("A",""), 
                    ggarrange( BC_WZA, BC_TC, ncol = 1, nrow = 2, legend = F, labels = c("B", "C")), common.legend = T),
                       BC_powerCurve, ncol = 2, widths = c(2,1), labels = c("","D"))

          



cline_snp_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/cline/sample40/3_0.003_3.12Loci.directionalSelection.sample.csv")
cline_WZA_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/cline/WZA_cline_sample40.csv")
cline_summary_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/cline/WZA_cline_sample40.summary.csv.summary.csv")


cline_allSNPs <- ggplot(data = cline_snp_file[cline_snp_file$maf > 0.05,], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  geom_point()+
  geom_point(data = cline_snp_file[(cline_snp_file$maf > 0.05)&(cline_snp_file$LA > 0.005),], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(italic(-log[10]*"(Kendall's "*tau*" p-value)")))+
  theme_bw()

cline_WZA <- ggplot(data = cline_WZA_file, aes(x = position/1e6, y = empR_Z, col = LA> 0.005))+
  geom_point()+
  geom_point(data = cline_WZA_file[(cline_WZA_file$LA > 0.005),], aes(x = position/1e6, y = empR_Z, col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(Z[W]))+  
  theme_bw()


cline_TC <- ggplot(data = cline_WZA_file, aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005))+
  geom_point()+
  geom_point(data = cline_WZA_file[(cline_WZA_file$LA > 0.005),], aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous("Top-candidate\nIndex")+
  theme_bw()


paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")


cline_powerCurve <- ggplot(data = cline_summary_file[cline_summary_file$rep == "mean",], aes(x = top, y = TP_WZA_empR))+
  geom_line(aes(col = "WZA"),lwd = 2)+
  geom_line(data = cline_summary_file[cline_summary_file$rep == "mean",], aes(x = top, y = TP_TC_uncor, col = "Top-Candidate"),lwd = 2)+
  geom_line(data = cline_summary_file[cline_summary_file$rep == "mean",], aes(x = top, y = TP_SNP_uncor, col = "SNP-based"),lwd = 2)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  scale_x_continuous("Number of Genes in Top Set")+
  scale_color_manual("Test", values = paletteCB)+
  theme_bw()

library(ggpubr)




cline_panel_plot <- ggarrange( ggarrange(cline_allSNPs, labels = c("E",""), 
                                      ggarrange( cline_WZA, cline_TC, ncol = 1, nrow = 2, legend = F, labels = c("F", "G")), common.legend = T),
                            cline_powerCurve, ncol = 2, widths = c(2,1), labels = c("","H"))








trunc_snp_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/trunc/sample40/10_0.003_2.12Loci.directionalSelection.sample.csv")
trunc_WZA_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/trunc/WZA_trunc_sample40.csv")
trunc_summary_file <- read.csv("~/work/GEA/simulations/directionalSelection/panelPlot/trunc/WZA_trunc_sample40.summary.csv.summary.csv")

trunc_snp_file[trunc_snp_file$maf > 0.05,]
trunc_allSNPs <- ggplot(data = trunc_snp_file[trunc_snp_file$maf > 0.05,], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  geom_point()+
  geom_point(data = trunc_snp_file[(trunc_snp_file$maf > 0.05)&(trunc_snp_file$LA > 0.005),], aes(x = position/1e6, y = -log10(geno_k_tau_p_value), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(italic(-log[10]*"(Kendall's "*tau*" p-value)")), limits = c(0,5))+
  theme_bw()

trunc_WZA <- ggplot(data = trunc_WZA_file, aes(x = position/1e6, y = empR_Z, col = LA> 0.005))+
  geom_point()+
  geom_point(data = trunc_WZA_file[(trunc_WZA_file$LA > 0.005),], aes(x = position/1e6, y = empR_Z, col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene",values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous(expression(Z[W]))+  
  theme_bw()


trunc_TC <- ggplot(data = trunc_WZA_file, aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005))+
  geom_point()+
  geom_point(data = trunc_WZA_file[(trunc_WZA_file$LA > 0.005),], aes(x = position/1e6, y = -log10(top_candidate_p), col = LA> 0.005 ))+
  scale_color_manual("Locally Adaptive Gene", values = c("darkgrey","black"))+
  scale_x_continuous("Position in Genome (Mbp)")+
  scale_y_continuous("Top-candidate\nIndex")+
  theme_bw()


paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")


trunc_powerCurve <- ggplot(data = trunc_summary_file[trunc_summary_file$rep == "mean",], aes(x = top, y = TP_WZA_empR))+
  geom_line(aes(col = "WZA"),lwd = 2)+
  geom_line(data = trunc_summary_file[trunc_summary_file$rep == "mean",], aes(x = top, y = TP_TC_uncor, col = "Top-Candidate"),lwd = 2)+
  geom_line(data = trunc_summary_file[trunc_summary_file$rep == "mean",], aes(x = top, y = TP_SNP_uncor, col = "SNP-based"),lwd = 2)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  scale_x_continuous("Number of Genes in Top Set")+
  scale_color_manual("Test", values = paletteCB)+
  theme_bw()

library(ggpubr)




trunc_panel_plot <- ggarrange( ggarrange(trunc_allSNPs, labels = c("I",""), 
                                         ggarrange( trunc_WZA, trunc_TC, ncol = 1, nrow = 2, legend = F, labels = c("J", "K")), common.legend = T),
                               trunc_powerCurve, ncol = 2, widths = c(2,1), labels = c("","L"))




png("~/work/GEA/simulations/Plots/ThreeRowDemoPlot.png", 
    units = "in", 
    width = 12, 
    height = 12,
    res = 300)
ggarrange( BC_panel_plot,
           cline_panel_plot,
           trunc_panel_plot,
           nrow = 3, ncol = 1)
dev.off()











 # Mike suggested removing the line plots from this figure, so here's a version wthout them

BC_panel_plot_noLine <-   annotate_figure( ggarrange(BC_allSNPs, BC_WZA, BC_TC, ncol = 3, nrow = 1, legend = "right", labels = c("I","J", "K"), common.legend = T), top = "BC Map")
cline_panel_plot_noLine <-   annotate_figure( ggarrange(cline_allSNPs, cline_WZA, cline_TC, ncol = 3, nrow = 1, legend = "right", labels = c("I","J", "K"), common.legend = T), top = "Gradient Map")
trunc_panel_plot_noLine <-   annotate_figure( ggarrange(trunc_allSNPs, trunc_WZA, trunc_TC, ncol = 3, nrow = 1, legend = "right", labels = c("I","J", "K"), common.legend = T), top = "Truncated Map")

png("~/work/GEA/simulations/Plots/ThreeRowDemoPlot_noLine.png", 
    units = "in", 
    width = 12, 
    height = 12,
    res = 300)
ggarrange( BC_panel_plot_noLine,
           cline_panel_plot_noLine,
           trunc_panel_plot_noLine,
           nrow = 3, ncol = 1)
dev.off()

