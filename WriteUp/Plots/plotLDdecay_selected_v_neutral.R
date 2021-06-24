rm(list= ls())

s_1 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_10_0.003_1.LD.csv", header = F)
s_1$rep <- 1
s_2 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_15_0.003_1.LD.csv", header = F)
s_2$rep <- 2
s_3 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_2_0.003_1.LD.csv", header = F)
s_3$rep <- 3
s_4 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_22_0.003_1.LD.csv", header = F)
s_4$rep <- 4
s_5 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_23_0.003_1.LD.csv", header = F)
s_5$rep <- 5
s_6 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_24_0.003_1.LD.csv", header = F)
s_6$rep <- 6
s_7 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_25_0.003_1.LD.csv", header = F)
s_7$rep <- 7
s_8 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_3_0.003_1.LD.csv", header = F)
s_8$rep <- 8
s_9 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_4_0.003_1.LD.csv", header = F)
s_9$rep <- 9
s_10 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/selected/selection_5_0.003_1.LD.csv", header = F)
s_10$rep <- 10

LD_s<-rbind( s_1,
             s_2,
             s_3,
             s_4,
             s_5,
             s_6,
             s_7,
             s_8,
             s_9,
             s_10)



n_1 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/10_0.003_1.LD.csv", header = F)
n_1$rep <- 1
n_2 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/15_0.003_1.LD.csv", header = F)
n_2$rep <- 2
n_3 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/2_0.003_1.LD.csv", header = F)
n_3$rep <- 3
n_4 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/22_0.003_1.LD.csv", header = F)
n_4$rep <- 4
n_5 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/23_0.003_1.LD.csv", header = F)
n_5$rep <- 5
n_6 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/24_0.003_1.LD.csv", header = F)
n_6$rep <- 6
n_7 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/25_0.003_1.LD.csv", header = F)
n_7$rep <- 7
n_8 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/3_0.003_1.LD.csv", header = F)
n_8$rep <- 8
n_9 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/4_0.003_1.LD.csv", header = F)
n_9$rep <- 9
n_10 <- read.csv("work/GEA/simulations/directionalSelection/G.2.3/LD_decay/neutral/5_0.003_1.LD.csv", header = F)
n_10$rep <- 10

library(ggplot2)
LD_n<-rbind( n_1,
             n_2,
             n_3,
             n_4,
             n_5,
             n_6,
             n_7,
             n_8,
             n_9,
             n_10)

LD_s$id <- paste(LD_s$V1, "_", LD_s$rep, sep ="")
LD_n$id <- paste(LD_n$V1, "_", LD_s$rep, sep ="")

ggplot()+
  
  geom_smooth(data = LD_s, aes(x = V2, y = V3), col = "red", fill = "red")+
  geom_smooth(data = LD_n, aes(x = V2, y = V3), col = "blue", fill = "blue")

