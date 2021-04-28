rm(list = ls())

library(ggplot2)

######## Plot effect size Distribution

bc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/BC_Map_sampled.WZA.csv")
cline_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/cline_sampled.WZA.csv")
trunc_10000_Vs192 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/Vs192/analysisFiles/trunc_sampled.WZA.csv")

bc_10000_Vs192$map <- "BC Map"
cline_10000_Vs192$map <- "Cline"
trunc_10000_Vs192$map <- "Truncated"
stablising_data <- rbind( bc_10000_Vs192, cline_10000_Vs192, trunc_10000_Vs192 )
stablising_data$model <- "Stabilising Selection"

bc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/BC_Map_sampled.WZA.csv")
cline_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/cline_sampled.WZA.csv")
trunc_10000_s0.003 <- read.csv("~/work/GEA/simulations/directionalSelection/G.2.3/BayPassResults/s0.003/analysisFiles/trunc_sampled.WZA.csv")


bc_10000_s0.003$map <- "BC Map"
cline_10000_s0.003$map <- "Cline"
trunc_10000_s0.003$map <- "Truncated"
directional_data <- rbind( bc_10000_s0.003, cline_10000_s0.003, trunc_10000_s0.003 )
directional_data$model <- "Directional Selection"

effect_size_data <- rbind( stablising_data, directional_data)


effectSizeDistribution_dir <- ggplot()+
  geom_histogram( data=directional_data[directional_data$LA > 0, ], aes( x= LA ), bins = 50, fill = "#91cf60", col = "black")+
  geom_vline(xintercept = 0.005, col = "red", lty = 2)+
  facet_wrap(model~ map, scales = "free")+
  xlab("Covariance of Phenotype and Environment")+
  ylab("Count")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )

effectSizeDistribution_stab <- ggplot()+
  geom_histogram( data=stablising_data[stablising_data$LA > 0, ], aes( x= LA ), bins = 60, fill = "#fc8d59", col = "black")+
  geom_vline(xintercept = 0.005, col = "red", lty = 2)+
  facet_wrap(model~ map, scales = "free")+
  xlab("Covariance of Phenotype and Environment")+
  ylab("Count")+
  theme_bw()+
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 12)
  )


pdf("~/work/GEA/simulations/Plots/effectSizeDistributionPlot.pdf", width = 12, height = 5)
print(effectSizeDistribution)  
dev.off()



BC_Map_phenotypes <- read.csv("/media/booker/HOWDY/G.2.3/Vs192/BC_Map/1_0.5_192.phen.txt")
BC_Map_phenotypes$map <-"BC Map"
cline_phenotypes <- read.csv("/media/booker/HOWDY/G.2.3/Vs192/cline/1_0.5_192.phen.txt")
cline_phenotypes$map <-"Gradient"
trunc_phenotypes <- read.csv("/media/booker/HOWDY/G.2.3/Vs192/trunc/1_0.5_192.phen.txt")
trunc_phenotypes$map <- "Truncated"

phen <- rbind(BC_Map_phenotypes, cline_phenotypes, trunc_phenotypes)

indPhens <- ggplot(data = phen , aes( x = opt, y = phen, fill = phen))+
  geom_jitter(shape = 21, height = 0, width = 0.2)+
  scale_fill_gradient2(low = "#CC6600", high = "#3399CC", mid = "white")+
  facet_grid(~map)+
  scale_x_continuous("Phenotypic Optimum")+
  scale_y_continuous("Phenotype", limits = c(-7,7))+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none"
    )

## Use the aggregation function to make a quick summary of the data.frame
temp <- aggregate(phen ~ map+pop+opt, phen, FUN = "mean")

phenMeans <- ggplot(data = temp , aes( x = opt, y = phen, fill =phen))+
  #stat_summary(fun.y = "mean", geom = "point", shape = 21)+
  geom_jitter(shape = 21, height = 0, width = 0.2)+
  scale_fill_gradient2(low = "#CC6600", high = "#3399CC", mid = "white")+
  facet_grid(~map)+
  scale_x_continuous("Phenotypic Optimum")+
  scale_y_continuous("Population Mean Phenotype", limits = c(-7,7))+
  theme_bw()+
  theme(
    panel.spacing.y = unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none"
  )

png("~/work/GEA/simulations/Plots/PhenotypePlot.png", units = "in", width = 10, height = 6, res = 300)
ggarrange(indPhens, phenMeans, nrow = 2, ncol = 1)
dev.off()
