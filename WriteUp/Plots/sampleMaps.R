rm(list = ls())

n40_samples <- c(119, 15, 148, 106, 104, 111, 18, 76, 50, 173, 166, 177, 101, 1, 153, 35, 149, 137, 40, 121, 110, 124, 89, 113, 98, 115, 51, 47, 96, 8, 120, 152, 91, 60, 31, 117, 143, 53, 49, 81)

bc <- read.csv("~/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt",header = F)

bc$lat <- rep(1:14,14) 
bc$lon <- rep(1:14,each = 14) 

library(ggplot2)

ggplot(data = bc, aes(x = lon, y = lat, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = bc[n40_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

bc_n20_samples <- c(153, 173, 143,  18,  53,  60, 152, 148,  15,  47, 166,  51, 121,  50, 149, 101, 111,  76,   1, 110)
ggplot(data = bc, aes(x = lon, y = lat, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = bc[bc_n20_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

bc_n10_samples <- c( 1,  51, 149, 110, 148,  15, 101, 153,  47, 166)
ggplot(data = bc, aes(x = lon, y = lat, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = bc[bc_n10_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()







## Gradient

cline <- read.csv("~/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt",header = F)

cline$lat <- rep(1:14,14) 
cline$lon <- rep(1:14,each = 14) 

library(ggplot2)

ggplot(data = cline, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = cline[n40_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

cline_n20_samples <- c(153, 173, 143,  18,  53,  60, 152, 148,  15,  47, 166,  51, 121,  50, 149, 101, 111,  76,   1, 110)
ggplot(data = cline, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = cline[cline_n20_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

cline_n10_samples <- c( 1,  51, 149, 110, 148,  15, 101, 153,  47, 166)
ggplot(data = cline, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = cline[cline_n10_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()




## Truncated

trunc <- read.csv("~/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt",header = F)

trunc$lat <- rep(1:14,14) 
trunc$lon <- rep(1:14,each = 14) 

library(ggplot2)

ggplot(data = trunc, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = trunc[n40_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

trunc_n20_samples <- c(153, 173, 143,  18,  53,  60, 152, 148,  15,  47, 166,  51, 121,  50, 149, 101, 111,  76,   1, 110)
ggplot(data = trunc, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = trunc[trunc_n20_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()

trunc_n10_samples <- c( 1,  51, 149, 110, 148,  15, 101, 153,  47, 166)
ggplot(data = trunc, aes(x = lat, y = lon, fill = V1))+
  geom_point(shape = 22, size = 10, col = "grey")+
  geom_point(data = trunc[trunc_n10_samples,], shape = 2, size = 3)+
  scale_fill_gradient2(high = "#3399CC", low = "#CC6600", mid = "white")+
  guides(col= F, fill = F)+
  theme_void()








