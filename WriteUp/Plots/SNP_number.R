rm(list = ls())
library(ggplot2)
library(ggpubr)

dir_trunc_2SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/trunc_sampled.WZA.n2.summary.csv")
dir_trunc_2SNPs$snps <- "WZA - 2 SNPs"

dir_trunc_4SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/trunc_sampled.WZA.n4.summary.csv")
dir_trunc_4SNPs$snps <- "WZA - 4 SNPs"

dir_trunc_8SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/trunc_sampled.WZA.n8.summary.csv")
dir_trunc_8SNPs$snps <- "WZA - 8 SNPs"

dir_trunc_w10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.summary.csv")
dir_trunc_w10000$snps <- "WZA - f10,000bp Window"

SNP_num_trunc <- rbind(dir_trunc_2SNPs, dir_trunc_4SNPs, dir_trunc_8SNPs)


SNP_plot_trunc <- ggplot(SNP_num_trunc[SNP_num_trunc$rep == "mean",], aes( x = top, y = TP_WZA_uncor, col = snps))+
  geom_line(lwd= 1.2)+
  geom_line(data=dir_trunc_w10000[dir_trunc_w10000$rep == "mean",], aes( x = top, y =TP_SNP_uncor, col = "Single SNP") , lwd= 1.2)+
  geom_line(data=dir_trunc_w10000[dir_trunc_w10000$rep == "mean",], aes( x = top, y =TP_WZA_uncor, col = "WZA - As in main text") , lwd= 1.2)+
  theme_bw()+
  ggtitle("Truncated")+
  scale_colour_brewer( palette = "Set1", guide = guide_legend(reverse = T))+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Number of Genes in Top Set")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    legend.title = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


dir_cline_2SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/cline_sampled.WZA.n2.summary.csv")
dir_cline_2SNPs$snps <- "WZA - 2 SNPs"

dir_cline_4SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/cline_sampled.WZA.n4.summary.csv")
dir_cline_4SNPs$snps <- "WZA - 4 SNPs"

dir_cline_8SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/cline_sampled.WZA.n8.summary.csv")
dir_cline_8SNPs$snps <- "WZA - 8 SNPs"

dir_cline_w10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled.WZA.summary.csv")
dir_cline_w10000$snps <- "WZA - 10,000bp Window"

SNP_num_cline <- rbind(dir_cline_2SNPs, dir_cline_4SNPs, dir_cline_8SNPs)

SNP_plot_cline <- ggplot(SNP_num_cline[SNP_num_cline$rep == "mean",], aes( x = top, y = TP_WZA_uncor, col = snps))+
  geom_line(lwd= 1.2)+
  geom_line(data=dir_cline_w10000[dir_cline_w10000$rep == "mean",], aes( x = top, y =TP_SNP_uncor, col = "Single SNP") , lwd= 1.2)+
  geom_line(data=dir_trunc_w10000[dir_trunc_w10000$rep == "mean",], aes( x = top, y =TP_WZA_uncor, col = "WZA - As in main text") , lwd= 1.2)+
  theme_bw()+
  ggtitle("Gradient")+
  scale_colour_brewer( palette = "Set1", guide = guide_legend(reverse = T))+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Number of Genes in Top Set")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    legend.title = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

dir_BC_Map_2SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/BC_Map_sampled.WZA.n2.summary.csv")
dir_BC_Map_2SNPs$snps <- "WZA - 2 SNPs"

dir_BC_Map_4SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/BC_Map_sampled.WZA.n4.summary.csv")
dir_BC_Map_4SNPs$snps <- "WZA - 4 SNPs"

dir_BC_Map_8SNPs <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/SNP_number_s0.003/BC_Map_sampled.WZA.n8.summary.csv")
dir_BC_Map_8SNPs$snps <- "WZA - 8 SNPs"

dir_BC_Map_w10000 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled.WZA.summary.csv")
dir_BC_Map_w10000$snps <- "WZA - 10,000bp Window"

SNP_num_BC_Map <- rbind(dir_BC_Map_2SNPs, dir_BC_Map_4SNPs, dir_BC_Map_8SNPs)

SNP_plot_BC_Map <- ggplot(SNP_num_BC_Map[SNP_num_BC_Map$rep == "mean",], aes( x = top, y = TP_WZA_uncor, col = snps))+
  geom_line(lwd= 1.2)+
  geom_line(data=dir_BC_Map_w10000[dir_BC_Map_w10000$rep == "mean",], aes( x = top, y =TP_SNP_uncor, col = "Single SNP") , lwd= 1.2)+
  geom_line(data=dir_trunc_w10000[dir_trunc_w10000$rep == "mean",], aes( x = top, y =TP_WZA_uncor, col = "WZA - As in main text") , lwd= 1.2)+
  theme_bw()+
  ggtitle("BC Map")+
  scale_colour_brewer( palette = "Set1", guide = guide_legend(reverse = T))+
  scale_y_continuous("Proportion of True Positives Detected", expand = c(0,0))+
  scale_x_continuous("Number of Genes in Top Set")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    legend.title = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )







pdf("~/work/GEA/simulations/Plots/SNP_number.pdf", height = 5, width = 14, onefile = F, labels = "AUTO")
ggarrange(SNP_plot_BC_Map,SNP_plot_cline, SNP_plot_trunc,nrow= 1, ncol = 3, common.legend = T, legend = "right") 
dev.off()


## Make a plot showing the SNPs in a genomic region with a strong outlier for each map

BC_Map_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.csv")
cline_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/cline_sampled.WZA.csv")
trunc_wza <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.csv")

BC_Map_MAX <- BC_Map_wza[order(BC_Map_wza$kendall_Z, decreasing = TRUE),]
head(BC_Map_MAX)
BC_Map_MAX$rep
BC_Map_MAX$gene

## Let's use Gene 726 from rep 4

BC_Map_instance<-read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled/4_0.003_directional_d40n50_i10000.csv")

BC_Map_rep <- BC_Map_wza[BC_Map_wza$rep== 4,]
TEMP <- BC_Map_instance[BC_Map_instance$gene %in% paste("gene",726, sep = ""),]
min(TEMP$position)
max(TEMP$position)

BC_Map_instance_SNPs <- ggplot(data = BC_Map_instance[BC_Map_instance$gene %in% paste("gene",(726-6):(726+6), sep = ""),], aes( x = position/1e6, y = -log10(geno_k_tau_p_value), col = pbar_qbar))+
  geom_point()+
  #  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((726-6):(726+7)), lty = 2, alpha = 0.1)+
  ggtitle("BC Map")+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  scale_y_continuous(expression(-log[10]*"("*italic("p-value")*")"),limits = c(0,11))+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


BC_Map_instance_WZA <- ggplot(data = BC_Map_rep[BC_Map_rep$gene %in% paste("gene",(726-6):(726+6), sep = ""),])+
  geom_point(aes(x = position/1e6, y = empR_Z), col = "plum4",size = 4)+
  #  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((726-6):(726+7)), lty = 2, alpha = 0.1)+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  scale_y_continuous(expression(italic(Z[W])), limits = c(-3,12))+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

BC_Map_plot <- ggarrange(
  ggarrange(BC_Map_instance_SNPs+theme(legend.position = "none"),BC_Map_instance_WZA,ncol = 1,nrow = 2,align = "v"),
  ggarrange(get_legend(BC_Map_instance_SNPs), ncol = 1,nrow = 2),
  widths = c(5,1)
)



cline_MAX <- cline_wza[order(cline_wza$kendall_Z, decreasing = TRUE),]
head(cline_MAX)
cline_MAX$rep
cline_MAX$gene

cline_rep <- cline_wza[cline_wza$rep== 7,]

cline_instance<-read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled/7_0.003_directional_d40n50_i10000.csv")

cline_instance_SNPs <- ggplot(data = cline_instance[cline_instance$gene %in% paste("gene",(330-6):(330+6), sep = ""),], aes( x = position/1e6, y = -log10(geno_k_tau_p_value), col = pbar_qbar))+
  geom_point()+
  #  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((330-6):(330+7)), lty = 2, alpha = 0.1)+
  ggtitle("Gradient")+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  scale_y_continuous(expression(-log[10]*"("*italic("p-value")*")"),limits = c(0,11))+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


cline_instance_WZA <- ggplot(data = cline_rep[cline_rep$gene %in% paste("gene",(330-6):(330+6), sep = ""),])+
  geom_point(aes(x = position/1e6, y = empR_Z), col = "plum4",size = 4)+
  #  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((330-6):(330+7)), lty = 2, alpha = 0.1)+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  scale_y_continuous(expression(italic(Z[W])), limits = c(-3,12))+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

cline_plot <- ggarrange(
  ggarrange(cline_instance_SNPs+theme(legend.position = "none"),cline_instance_WZA,ncol = 1,nrow = 2,align = "v"),
  ggarrange(get_legend(cline_instance_SNPs), ncol = 1,nrow = 2),
  widths = c(5,1)
)




trunc_MAX <- trunc_wza[order(trunc_wza$kendall_Z, decreasing = TRUE),]
head(trunc_MAX)
trunc_MAX$rep
trunc_MAX$gene

trunc_rep <- trunc_wza[trunc_wza$rep== 11,]

trunc_instance<-read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled/11_0.003_directional_d40n50_i10000.csv")
trunc_instance[trunc_instance$gene=="gene594",]
trunc_instance_SNPs <- ggplot(data = trunc_instance[trunc_instance$gene %in% paste("gene",(594-6):(594+6), sep = ""),], aes( x = position/1e6, y = -log10(geno_k_tau_p_value), col = pbar_qbar))+
  geom_point()+
#  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((594-6):(594+7)), lty = 2, alpha = 0.1)+
  ggtitle("Truncated")+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  scale_y_continuous(expression(-log[10]*"("*italic("p-value")*")"),limits = c(0,11))+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )


trunc_instance_WZA <- ggplot(data = trunc_rep[trunc_rep$gene %in% paste("gene",(594-6):(594+6), sep = ""),])+
  geom_point(aes(x = position/1e6, y = empR_Z), col = "plum4",size = 4)+
  #  geom_smooth(col = "red",aes(weight = pbar_qbar), span=0.5)+
  geom_vline(xintercept = 0.01*((594-6):(594+7)), lty = 2, alpha = 0.1)+
  scale_x_continuous("Position in Chromosome (Mbp)")+
  scale_y_continuous(expression(italic(Z[W])), limits = c(-3,12))+
  scale_color_gradient(expression(italic(bar(p)*" "*bar(q))),low = "grey", high = "black")+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

trunc_plot <- ggarrange(
  ggarrange(trunc_instance_SNPs+theme(legend.position = "none"),trunc_instance_WZA,ncol = 1,nrow = 2,align = "v"),
  ggarrange(get_legend(trunc_instance_SNPs), ncol = 1,nrow = 2),
  widths = c(5,1)
)


pdf("~/work/GEA/simulations/Plots/cline_plot_demo.pdf", height = 7, width = 10, onefile = F)
print(cline_plot)
dev.off()

pdf("~/work/GEA/simulations/Plots/BC_Map_plot_demo.pdf", height = 7, width = 10, onefile = F)
print(BC_Map_plot)
dev.off()

pdf("~/work/GEA/simulations/Plots/trunc_plot_demo.pdf", height = 7, width = 10, onefile = F)
print(trunc_plot)
dev.off()

pdf("~/work/GEA/simulations/Plots/all_maps_plot_demo.pdf", height = 21, width = 10, onefile = F)
ggarrange(BC_Map_plot, cline_plot, trunc_plot, nrow = 3, ncol = 1, labels = "AUTO")
dev.off()
