rm(list = ls ())


bc_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_n10.WZA.csv.summary.csv")
bc_n10$source <- "n = 10"

bc_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_n20.WZA.summary.csv")
bc_n20$source <- "n = 20"

bc_n40 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled.WZA.summary.csv")
bc_n40$source <- "n = 40"

bc_all <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map.WZA2.summary.csv")
bc_all$source <- "Entire Population"

library(ggplot2)

bc <- rbind( bc_all, bc_n10, bc_n20, bc_n40)
bc$map <- "BC Map"


cline_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled_n10.WZA.summary.csv")
cline_n10$source <- "n = 10"

cline_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled_n20.WZA.summary.csv")
cline_n20$source <- "n = 20"

cline_n40 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_sampled.WZA.summary.csv")
cline_n40$source <- "n = 40"

cline_all <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline.WZA.summary.csv")
cline_all$source <- "Entire Population"

library(ggplot2)

cline <- rbind( cline_all, cline_n10, cline_n20, cline_n40)
cline$map <- "Gradient"


trunc_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled_n10.WZA.summary.csv")
trunc_n10$source <- "n = 10"

trunc_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled_n20.WZA.summary.csv")
trunc_n20$source <- "n = 20"

trunc_n40 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc_sampled.WZA.summary.csv")
trunc_n40$source <- "n = 40"

trunc_all <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/trunc.WZA.summary.csv")
trunc_all$source <- "Entire Population"

library(ggplot2)

trunc <- rbind( trunc_all, trunc_n10, trunc_n20, trunc_n40)
trunc$map <- "Truncated"

gea <- rbind( bc, cline, trunc)


paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")


sampleSizeComparison <- ggplot(data = gea[gea$rep == "mean",])+
  geom_line(aes(x = top, y = TP_SNP_uncor, col = "SNP-based"), lwd = 1)+
  geom_line(aes(x = top, y = TP_TC_uncor, col = "Top-Candidate"), lwd = 1)+
  geom_line(aes(x = top, y = TP_WZA_empR, col = "WZA"), lwd = 1)+
  facet_grid(source~map)+
  scale_color_manual("Test", values = paletteCB)+
  scale_y_continuous("Proportion of True Positives Detected", limits = c(0,1))+
  theme_bw()+
  scale_x_continuous("Top # Genes")+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

pdf("~/work/GEA/simulations/Plots/sampleSizeComparison.pdf", width = 8, height = 10)
print(sampleSizeComparison)
dev.off()
