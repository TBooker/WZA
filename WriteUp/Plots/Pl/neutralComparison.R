rm(list = ls())

library(ggplot2)
library(ggpubr)

bcMap <- read.csv("/media/booker/HOWDY/GEA/G.0.9/BCmap_neutral.WZA.csv")

BC<-ggplot(data = bcMap, aes(x = SNPs, y = Z))+
  geom_point(width = 0.1, alpha=0.4)+
  theme_bw()+
  ggtitle("Stepping-Stone Model (BC Map)")+
  scale_y_continuous(limits = c(-5,15))

BC_QQ <- ggplot(data = bcMap, aes(sample = Z))+
  stat_qq( alpha=0.4)+
  stat_qq_line( alpha=0.4)+
  theme_bw()+
  scale_y_continuous("Sample",limits = c(-10,20))+
  scale_x_continuous("Theoretical",limits = c(-5,5))

BC_hist<- ggplot(data = bcMap, aes(x = Z))+
  geom_histogram( alpha=0.75,col = "black", bins = 50)+
  theme_bw()+
  scale_x_continuous(limits = c(-10,20))+
  scale_y_continuous("Count")




islandMap <- read.csv("/media/booker/HOWDY/GEA/G.0.9/IslandModel_neutral.WZA.csv")

island<-ggplot(data = islandMap, aes(x = SNPs, y = Z))+
  geom_point(width = 0.1, alpha=0.4)+
  theme_bw()+
  ggtitle("Island Model")+
  scale_y_continuous(limits = c(-5,15))


island_QQ <- ggplot(data = islandMap, aes(sample = Z))+
  stat_qq( alpha=0.4)+
  stat_qq_line( alpha=0.4)+
  theme_bw()+
  scale_y_continuous("Sample",limits = c(-10,20))+
  scale_x_continuous("Theoretical",limits = c(-5,5))

island_hist<-ggplot(data = islandMap, aes(x = Z))+
    geom_histogram( alpha=0.75,col = "black", bins = 50)+
    theme_bw()+
    scale_x_continuous(limits = c(-10,20))+
    scale_y_continuous("Count")

  
  
  
ggarrange(island, BC, island_QQ, BC_QQ, island_hist, BC_hist, nrow = 3, ncol = 2)
