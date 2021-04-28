library(dplyr)

setwd("~/UBC/GEA/WZA/SimulationStudy/directionalSelection/G.2.23/isolation_by_distance/")
temp = list.files(pattern = "*trees")
myfiles = lapply(temp, read.csv)

FstDF <- bind_rows(myfiles, .id = "column_label")
FstDF$source <- NA
FstDF[FstDF$dist ==-1,]$source = "Population Mean"
FstDF[FstDF$dist >0,]$source = "Isolation-by-Distance"

library(ggplot2)
meanFst = mean(FstDF[FstDF$dist == -1,]$Fst )

FstPlot <- ggplot(data = FstDF[FstDF$source == "Isolation-by-Distance",] , aes(x = dist, y = Fst))+
  geom_hline(aes( yintercept = meanFst) , lty = 2, col = "red")+
  geom_smooth()+
  geom_point()+
  scale_x_continuous("Distance Between Demes", breaks = 1:14)+
  scale_y_continuous(expression(F[ST]), limits = c(0,0.15))+
  theme_bw()


setwd("~/UBC/GEA/WZA/SimulationStudy/directionalSelection/G.2.23/LD_decay/")
LD_temp = list.files(pattern = "*csv")
myfiles_LD = lapply(LD_temp, read.csv, header = F)

LDDF <- bind_rows(myfiles_LD, .id = "column_label")


LD_Plot <- ggplot(data = LDDF , aes(x = V1, y = V2, group = column_label))+
  geom_smooth()+
  scale_y_continuous(expression(r^2))+
  scale_x_continuous("Distance between sites (bp)")+
  theme_bw()

library(ggpubr)

png("~/UBC/GEA/WZA/SimulationStudy/directionalSelection/SummaryStats.png", width = 8, height = 5, units = "in", res = 300)
ggarrange( FstPlot, LD_Plot, labels = "AUTO", nrow = 1, ncol = 2)
dev.off()
