rm(list = ls())

## Let's make a plot of WZA results from Lodgepole Pine DD0 data

pine <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine_BF.csv")

library(ggplot2)

hits_to_highlight <- pine[(-log10(pine$top_candidate_p) < 10)&(pine$Z>10),]$contig

max_Z <- pine[pine$Z == max(pine$Z),]$contig

### Read in the raw p-values for four contigs
#raw <- read.csv("~/UBC/Convergence/Yeaman_Data/doi_10.5061_dryad.0t407__v1/yeaman_pine_spruce_convergence_archive2/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_raw", sep = "\t")
#scaffold750644_27863_31476
raw_1 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/tscaffold2670_552417_570968.tsv", header = T, comment.char = "&")
Z_1 <- pine[pine$contig == "tscaffold2670_552417_570968", ]
#tscaffold1962_233534_237733
raw_2 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/tscaffold1962_233534_237733.tsv", header = T, comment.char = "&")
Z_2 <- pine[pine$contig == "tscaffold1962_233534_237733", ]
#scaffold878015.2_134352_141334
raw_3 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/scaffold878015.2_134352_141334.tsv", header = T, comment.char = "&")
Z_3 <- pine[pine$contig == "scaffold878015.2_134352_141334", ]
#tscaffold3400_439096_449060
raw_4 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/tscaffold3400_439096_449060.tsv", header = T, comment.char = "&")
Z_4 <- pine[pine$contig == "tscaffold3400_439096_449060", ]

x = 10

colour_palette <- c("#e66101",  "#fdb863",  "#b2abd2",  "#5e3c99")
colour_palette <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')

corre <-cor( -log10(pine$top_candidate_p), y = pine$Z, method = "spearman")

cor_text = bquote("Spearman's" ~ rho ~ "=" ~ .(corre))

png("~/UBC/GEA/WZA/dataAnalsis/topCandidate_versus_WZA_DD_0.png", units = "in", height = 6, width = 6, res = 300)
ggplot( data = pine, aes(x = -log10(top_candidate_p), y = Z))+
  geom_point(alpha = 0.5)+
  scale_y_continuous(expression(Z[W]))+
  scale_x_continuous(expression(-log[10]*"(Top-Candidate Test "*italic("p")*"-value)"))+
  annotate("text", label = cor_text, x = 60, y = 0)+
  theme_bw()
dev.off()

part_1 <- ggplot( data = pine, aes(x = -log10(top_candidate_p), y = Z))+
  geom_point(alpha = 0.5)+
  geom_point( data = Z_1, aes(x = -log10(top_candidate_p), y = Z), fill = colour_palette[1], size = 3, shape = 21)+
  geom_point( data = Z_2, aes(x = -log10(top_candidate_p), y = Z), fill = colour_palette[2], size = 3, shape = 21)+
  geom_point( data = Z_3, aes(x = -log10(top_candidate_p), y = Z), fill = colour_palette[3], size = 3, shape = 21)+
  geom_point( data = Z_4, aes(x = -log10(top_candidate_p), y = Z), fill = colour_palette[4], size = 3, shape = 21)+
  scale_y_continuous(expression(Z[W]))+
  scale_x_continuous(expression(-log[10]*"(Top-Candidate Test "*italic("p")*"-value)"))+
  theme_bw()

MAF_plot_1 <- ggplot( data = raw_1, aes(x = minor_freq, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[1],  shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_1 <- ggplot( data = raw_1, aes(x = pos_gcontig, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[1],  shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Position in Contig")


MAF_plot_2 <- ggplot( data = raw_2, aes(x = minor_freq, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[2], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_2 <- ggplot( data = raw_2, aes(x = pos_gcontig, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[2], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173)), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Position in Contig")


MAF_plot_3 <- ggplot( data = raw_3, aes(x = minor_freq, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[3], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Minor Allele Frequency")


POS_plot_3 <- ggplot( data = raw_3, aes(x = pos_gcontig, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[3], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Position in Contig")


MAF_plot_4 <- ggplot( data = raw_4, aes(x = minor_freq, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[4], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_4 <- ggplot( data = raw_4, aes(x = pos_gcontig, y = log10(DD_0)))+
  geom_point(alpha = 0.9, fill = colour_palette[4], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = log10(3.473173) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(log[10]*"(Bayes Factor)"),limits = c(-1,30))+
  scale_x_continuous("Position in Contig")

library(ggpubr)


MAF_part_2 <- ggarrange( MAF_plot_1, MAF_plot_2, MAF_plot_3, MAF_plot_4, nrow = 2, ncol = 2, labels = c("B","C","D","E"))
POS_part_2 <- ggarrange( POS_plot_1, POS_plot_2, POS_plot_3, POS_plot_4, nrow = 2, ncol = 2, labels = c("B","C","D","E"))

png("~/UBC/GEA/WZA/dataAnalsis/Z_v_MAF_DD0.png", res = 300, height = 10, width = 5, units= "in")
ggarrange(part_1, MAF_part_2, labels = c("A", ""), nrow = 2, heights = c(2,3))
dev.off()

png("~/UBC/GEA/WZA/dataAnalsis/Z_v_POS_DD0.png", res = 300, height = 10, width = 5, units= "in")
ggarrange(part_1, POS_part_2, labels = c("A", ""), nrow = 2, heights = c(2,3))
dev.off()



pine_bay <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine_BF.csv")
pine_raw <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine.csv")

pine_bay <- pine_bay[sort(pine_bay$contig),]
pine_raw <- pine_raw[sort(pine_raw$contig),]

both <- merge(pine_bay, pine_raw, by.x="contig", by.y = "gene")

par(mfrow = c(2,2))

plot(-log10(both$top_candidate_p.y), -log10(both$top_candidate_p.x),
     xlab = "Top-Candidate Index - Structure Corrected",
     ylab = "Top-Candidate Index - Uncorrected",
     pch = 1,
     cex = 0.5)

plot(-log10(both$top_candidate_p.y), -log10(both$top_candidate_p.x),
     xlab = "Top-Candidate Index - Structure Corrected",
     ylab = "Top-Candidate Index - Uncorrected",
     pch = 1,
     cex = 0.5,
     xlim = c(0,10),
     ylim = c(0,10))

plot(both$Z.y, both$Z.x,
     xlab = "Z - Structure Corrected",
     ylab = "Z - Uncorrected",
     pch = 1,
     cex = 0.5)

plot(both$Z.y, both$Z.x,
     xlab = "Z - Structure Corrected",
     ylab = "Z - Uncorrected",
     pch = 1,
     cex = 0.5,
     xlim = c(-10,10),
     ylim = c(-10,10))


