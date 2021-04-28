rm( list = ls() )

ld1 <- read.csv("work/GEA/simulations/neutralSummaryStats/1.200102.LD.csv", sep = ',', head = 0)
ld1$rep = 1
ld2 <- read.csv("work/GEA/simulations/neutralSummaryStats/10.170102.LD.csv", sep = ',', head = 0)
ld2$rep = 2
ld3 <- read.csv("work/GEA/simulations/neutralSummaryStats/16.170102.LD.csv", sep = ',', head = 0)
ld3$rep = 3
ld4 <- read.csv("work/GEA/simulations/neutralSummaryStats/17.170102.LD.csv", sep = ',', head = 0)
ld4$rep = 4
ld5 <- read.csv("work/GEA/simulations/neutralSummaryStats/19.170102.LD.csv", sep = ',', head = 0)
ld5$rep = 5
ld6 <- read.csv("work/GEA/simulations/neutralSummaryStats/21.170102.LD.csv", sep = ',', head = 0)
ld6$rep = 6
ld7 <- read.csv("work/GEA/simulations/neutralSummaryStats/22.170102.LD.csv", sep = ',', head = 0)
ld7$rep = 7
ld8 <- read.csv("work/GEA/simulations/neutralSummaryStats/2.170102.LD.csv", sep = ',', head = 0)
ld8$rep = 8
ld9 <- read.csv("work/GEA/simulations/neutralSummaryStats/20.170102.LD.csv", sep = ',', head = 0)
ld9$rep = 9
ld10 <- read.csv("work/GEA/simulations/neutralSummaryStats/25.170102.LD.csv", sep = ',', head = 0)
ld10$rep = 10

library(ggplot2)
str(ld1)
theory = function(Ne, r){
  return ( 1/(1 + 4*Ne * r) )
}

r <- 1e-7*((1:100)/10)
expectation <- data.frame(rec = r , dist = (1:100)/10)
expectation$r2<-theory(100*196, r*1000)
str(expectation)
ggplot(data = expectation, aes(x = dist, y = r2))+
  geom_line()

#ggplot(data = rbind(ld1,ld2,ld3,ld4,ld5,ld6,ld7,ld8,ld9,ld10), aes(x = V1/1000, y = V2, group = rep))+
ggplot()+
  #  geom_smooth(data = rbind(ld1,ld2,ld3,ld4,ld5,ld6,ld7,ld8,ld9,ld10), aes(x = V1/1000, y = V2, group = rep), span = 0.1, col = "black")+
    geom_smooth(span = 5, method = "loess",data = ld1, aes(x = V1/1000, y = V2, group = rep), col = "black")+
  #  geom_smooth(span = 0.01, col = "black")+
  geom_line(data = expectation, aes(x = dist, y = r2), col ="red")+
    scale_y_continuous(expression(italic(r^2)))+
  scale_x_continuous("Distance between SNPs (Kbp)")+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black"),
    plot.title  =  element_text(size = 12, angle= 0, colour = "black", hjust = 0.5)
  )
?geom_smooth
