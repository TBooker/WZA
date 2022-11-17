rm(list = ls())

bc <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/s0.003/BC_Map_WZA_baypass.csv")
bc$map <- "BC Map"
cline <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/s0.003/cline_WZA_baypass.csv")
cline$map <- "Cline"

wza <- rbind(bc, cline)

library(ggplot2)
library(ggpubr)

bc_density_bayPass <- ggplot(bc, aes( x = Z_bayPass_, fill = LA > 0) )+
  geom_density(alpha= 0.3)+
  theme_bw()

bc_density_kendall <- ggplot(bc, aes( x = Z_kendall_, fill = LA > 0) )+
  geom_density(alpha= 0.3)+
  theme_bw()

bc_scatter <- ggplot(bc, aes( x = Z_kendall_, y = Z_bayPass_, col = LA > 0) )+
  geom_point(alpha= 0.4)+
  geom_density_2d( col = "black")+
  theme_bw()

cline_density_bayPass <- ggplot(cline, aes( x = Z_bayPass_, fill = LA > 0) )+
  geom_density(alpha= 0.3)+
  theme_bw()

cline_density_kendall <- ggplot(cline, aes( x = Z_kendall_, fill = LA > 0) )+
  geom_density(alpha= 0.3)+
  theme_bw()

cline_scatter <- ggplot(cline, aes( x = Z_kendall_, y = Z_bayPass_, col = LA > 0) )+
  geom_point(alpha= 0.4)+
  geom_density_2d( col = "black")+
  theme_bw()

ggarrange(bc_scatter, cline_scatter, bc_density_kendall,  cline_density_kendall , bc_density_bayPass, cline_density_bayPass,   nrow = 3, ncol = 2, common.legend = T)


cline_scatter <- ggplot(cline, aes( x = Z_bayPass_, y = -log10(BP_top_candidate_p), col = LA > 0) )+
  geom_point()+
  theme_bw()

cline_scatter <- ggplot(cline, aes( x = Z_kendall_, y = -log10(top_candidate_p), col = LA > 0) )+
  geom_point()+
  theme_bw()


ggplot(bc, aes( x = SNPs, y = hits, col = LA > 0) )+
  geom_point()+
  theme_bw()

str(bc)

ggplot(bc, aes( x = SNPs, y = Z_kendall_, col = LA > 0) )+
  geom_point()+
  theme_bw()
ggplot(bc, aes( x = SNPs, y = -log10(top_candidate_p), col = LA > 0) )+
  geom_point()+
  theme_bw()
