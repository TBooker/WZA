rm(list = ls())

## Let's make a plot of WZA results from Lodgepole Pine DD0 data

pine <- read.csv("~/UBC/GEA/WZA/dataAnalsis/DD_0_YeamanData_WZA_TC_empiricalP.csv")

mean(pine$Z)
sd(pine$Z)


cor.test(pine$Z, -log10(pine$top_candidate_p), method = "kendall")

library(ggplot2)
## Hits that are specific to WZA- grab 1 of these
WZA_hits_to_highlight <- pine[(-log10(pine$top_candidate_p) < 3)&(pine$Z>10),]$contig
WZA_hits_to_highlight[1]
## Hits that are specific to TC- grab 1 of these
TC_hits_to_highlight <- pine[(-log10(pine$top_candidate_p) > 6)&(pine$Z<5),]$contig
TC_hits_to_highlight[1]
## Hits that are found using either method - grab 2 of these
BOTH_hits_to_highlight <- pine[(-log10(pine$top_candidate_p) > 10)&(pine$Z>10),]$contig
as.character(BOTH_hits_to_highlight[1:2])

ScaffoldsToFocusOn <- c(as.character(WZA_hits_to_highlight[1]),
 as.character(TC_hits_to_highlight[1]),
 as.character(BOTH_hits_to_highlight[1:2]))

writeLines(ScaffoldsToFocusOn, "~/UBC/GEA/WZA/dataAnalsis/ScaffoldsToFocusOn.txt")

### Read in the raw p-values for the four contigs that are being focussed on
#raw <- read.csv("~/UBC/Convergence/Yeaman_Data/doi_10.5061_dryad.0t407__v1/yeaman_pine_spruce_convergence_archive2/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_raw", sep = "\t")
#scaffold750644_27863_31476
raw_1 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/scaffold133090_50563_54310.tsv", header = T, comment.char = "&")
Z_1 <- pine[pine$contig == "scaffold133090_50563_54310", ]
#tscaffold1962_233534_237733
raw_2 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/scaffold407823_54945_61669.tsv", header = T, comment.char = "&")
Z_2 <- pine[pine$contig == "scaffold407823_54945_61669", ]
#scaffold878015.2_134352_141334
raw_3 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/C32549516_106128_111892.tsv", header = T, comment.char = "&")
Z_3 <- pine[pine$contig == "C32549516_106128_111892", ]
#tscaffold3400_439096_449060
raw_4 <- read.table ("~/UBC/GEA/WZA/dataAnalsis/C32557620_55842_67865.tsv", header = T, comment.char = "&")
Z_4 <- pine[pine$contig == "C32557620_55842_67865", ]

x = 10

##Make a PNW appropriate colour pallete for the plot
colour_palette <- pnw_palette(name="Sunset2",n=4,type="discrete")

corre <-cor( -log10(pine$top_candidate_p), y = pine$Z, method = "spearman")

cor_text = bquote("Spearman's" ~ rho ~ "=" ~ .(corre))

png("~/UBC/GEA/WZA/dataAnalsis/topCandidate_versus_WZA_DD_0.png", units = "in", height = 6, width = 6, res = 300)
ggplot( data = pine, aes(x = -log10(top_candidate_p), y = Z))+
  geom_point(alpha = 0.5)+
  scale_y_continuous(expression(Z[W]))+
  scale_x_continuous(expression(-log[10]*"(Top-Candidate Test "*italic("p")*"-value)"))+
#  annotate("text", label = cor_text, x = 60, y = 0)+
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

MAF_plot_1 <- ggplot( data = raw_1, aes(x = minor_freq, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[1],  shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_1 <- ggplot( data = raw_1, aes(x = pos_gcontig, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[1],  shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Position in Contig")


MAF_plot_2 <- ggplot( data = raw_2, aes(x = minor_freq, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[2], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_2 <- ggplot( data = raw_2, aes(x = pos_gcontig, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[2], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Position in Contig")


MAF_plot_3 <- ggplot( data = raw_3, aes(x = minor_freq, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[3], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Minor Allele Frequency")


POS_plot_3 <- ggplot( data = raw_3, aes(x = pos_gcontig, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[3], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Position in Contig")


MAF_plot_4 <- ggplot( data = raw_4, aes(x = minor_freq, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[4], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
  scale_x_continuous("Minor Allele Frequency")

POS_plot_4 <- ggplot( data = raw_4, aes(x = pos_gcontig, y = -log10(DD_0_ep)))+
  geom_point(alpha = 0.9, fill = colour_palette[4], shape = 21)+
  theme_bw()+
  geom_hline(aes( yintercept = -log10(0.01) ), lty = 2, alpha = 0.6)+
  scale_y_continuous(expression(-log[10]*"(empirical p-value)"),limits = c(-1, 6))+
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


# 
# pine_bay <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine_BF.csv")
# pine_raw <- read.csv("~/UBC/GEA/WZA/dataAnalsis/AdapTreeData/DD_0_Pine.csv")
# 
# pine_bay <- pine_bay[sort(pine_bay$contig),]
# pine_raw <- pine_raw[sort(pine_raw$contig),]
# 
# both <- merge(pine_bay, pine_raw, by.x="contig", by.y = "gene")
# 
# par(mfrow = c(2,2))
# 
# plot(-log10(both$top_candidate_p.y), -log10(both$top_candidate_p.x),
#      xlab = "Top-Candidate Index - Structure Corrected",
#      ylab = "Top-Candidate Index - Uncorrected",
#      pch = 1,
#      cex = 0.5)
# 
# plot(-log10(both$top_candidate_p.y), -log10(both$top_candidate_p.x),
#      xlab = "Top-Candidate Index - Structure Corrected",
#      ylab = "Top-Candidate Index - Uncorrected",
#      pch = 1,
#      cex = 0.5,
#      xlim = c(0,10),
#      ylim = c(0,10))
# 
# plot(both$Z.y, both$Z.x,
#      xlab = "Z - Structure Corrected",
#      ylab = "Z - Uncorrected",
#      pch = 1,
#      cex = 0.5)
# 
# plot(both$Z.y, both$Z.x,
#      xlab = "Z - Structure Corrected",
#      ylab = "Z - Uncorrected",
#      pch = 1,
#      cex = 0.5,
#      xlim = c(-10,10),
#      ylim = c(-10,10))
# 

NormalFit <-  data.frame( dens =  dnorm( seq(-1000,1000)/100 , mean(pine$Z), sd(pine$Z) ),
                          Zw = seq(-1000,1000)/100)

Pine_Z_dist <- ggplot(data = pine, aes(x = Z))+
  geom_histogram(  aes(y=..density..),fill = colour_palette[2], bins = 80, col = "black")+
  geom_line(data = NormalFit, aes(x = Zw, y = dens))+
  scale_x_continuous(expression(italic(Z[W])))+
  scale_y_continuous(expression(italic("Density")))+
  theme_bw()

library(PNWColors)
pal <- pnw_palette(name="Sunset2",n=4,type="discrete")

pine_env <- read.csv("~/UBC/GEA/WZA/Yeaman_Data/Pine.Assn2.pops-envi.txt", sep = "\t")

DD0_dist <- ggplot(data = pine_env, aes( x  = DD_0))+
  geom_histogram(aes(y=..density..), fill = colour_palette[3], colour  = "black", bins = 20)+
  scale_y_continuous(expression(italic("Density")))+
  scale_x_continuous(expression(italic("Degree Days < 0")))+
  theme_bw()

pdf("~/UBC/GEA/WZA/dataAnalsis/pineAnalysis_SuppFig.pdf", height = 10, width = 6)
ggarrange(DD0_dist, Pine_Z_dist, ncol = 1, nrow = 2, labels = "AUTO")
dev.off()


x<-rnorm(100)
y <- x + rnorm(100)
cor(2*x,5*y)
