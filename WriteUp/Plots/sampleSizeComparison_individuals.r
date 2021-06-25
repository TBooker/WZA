rm(list = ls())

bc_d10_n5 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d10_n5.WZA.summary.csv")
bc_d10_n5$inds <- "5 Individuals Sampled"
bc_d10_n5$demes <- "10 Demes Sampled"

bc_d10_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d10_n10.WZA.summary.csv")
bc_d10_n10$inds <- "10 Individuals Sampled"
bc_d10_n10$demes <- "10 Demes Sampled"

bc_d10_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d10_n20.WZA.summary.csv")
bc_d10_n20$inds <- "20 Individuals Sampled"
bc_d10_n20$demes <- "10 Demes Sampled"

bc_d20_n5 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d20_n5.WZA.summary.csv")
bc_d20_n5$inds <- "5 Individuals Sampled"
bc_d20_n5$demes <- "20 Demes Sampled"

bc_d20_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d20_n10.WZA.summary.csv")
bc_d20_n10$inds <- "10 Individuals Sampled"
bc_d20_n10$demes <- "20 Demes Sampled"

bc_d20_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d20_n20.WZA.summary.csv")
bc_d20_n20$inds <- "20 Individuals Sampled"
bc_d20_n20$demes <- "20 Demes Sampled"

bc_d40_n5 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d40_n5.WZA.summary.csv")
bc_d40_n5$inds <- "5 Individuals Sampled"
bc_d40_n5$demes <- "40 Demes Sampled"

bc_d40_n10 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d40_n10.WZA.summary.csv")
bc_d40_n10$inds <- "10 Individuals Sampled"
bc_d40_n10$demes <- "40 Demes Sampled"

bc_d40_n20 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_sampled_d40_n20.WZA.summary.csv")
bc_d40_n20$inds <- "20 Individuals Sampled"
bc_d40_n20$demes <- "40 Demes Sampled"

library(ggplot2)

bc <- rbind( bc_d10_n5,
             bc_d10_n10,
             bc_d10_n20,
             bc_d20_n5,
             bc_d20_n10,
             bc_d20_n20,
             bc_d40_n5,
             bc_d40_n10,
             bc_d40_n20)
bc$map <- "BC Map"

library(ggplot2)

paletteCB <- c("#999999",
               "#56B4E9",
               "#D55E00")

bc$inds <- factor( bc$inds, 
                     levels = c( "5 Individuals Sampled",  "10 Individuals Sampled",  "20 Individuals Sampled"))

sampleSizeComparison <- ggplot(data = bc[bc$rep == "mean",])+
  geom_line(aes(x = top/1000, y = TP_SNP_uncor, col = "Kendall's tau"), lwd = 1)+
  geom_line(aes(x = top/1000, y = TP_TC_uncor, col = "Top-Candidate"), lwd = 1)+
  geom_line(aes(x = top/1000, y = TP_WZA_empR, col = "WZA"), lwd = 1)+
  facet_grid(demes~inds)+
  scale_color_manual("Method", values = paletteCB, labels = c(expression("WZA"*tau),
                                                              expression("Top-candidate"),
                                                              expression("Kendall's "* tau)),
                     breaks = c("WZA","Top-Candidate","Kendall's tau"))+
  scale_y_continuous("Proportion of true positives detected", limits = c(0,1))+
  theme_bw()+
  scale_x_continuous('Fraction of all genes in the top set')+
  coord_cartesian(ylim = c(0,1))+
  theme(
    panel.spacing.y = unit(1,"lines"),
    strip.background = element_blank(),
    #    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 30),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.text = element_text(hjust = 0)
  )

pdf( "~/work/GEA/simulations/Plots/sampleSizeComparison_individuals.pdf", width = 8, height = 6.5)
print(sampleSizeComparison)
dev.off()
