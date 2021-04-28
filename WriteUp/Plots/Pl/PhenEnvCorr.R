rm(list = ls())

corrs = list()
for (i in seq(61, 80)){
  temp = read.csv(paste("work/GEA/simulations/G.0.2/30101.",i,".checkMeOutChump.phen.txt", sep = ''))
  x = cor(temp$opt, temp$phen, method = "kendall")
  corrs[i - 40] = x
}

temp = read.csv(paste("~/work/GEA/simulations/G.0.2/30101.",i,".checkMeOutChump.phen.txt", sep = ''))

Correlations = data.frame(Correlation = unlist( corrs ), rep = seq(41, 60), parameters = "Vs=100")

library(ggplot2)
library(ggpubr)

phen <- read.csv('~/work/GEA/simulations/G.0.2/30101.61.checkMeOutChump.phen.txt')
cor.test(phen$phen, phen$opt, method = "kendall")

dataPlot <- ggplot(data = phen, aes( x = as.factor(opt), y = phen, fill = phen))+
  geom_jitter( shape = 21, col= 'black')+
  scale_y_continuous("Phenotype", limits = c(-10,10))+
  scale_x_discrete("Optimum")+
  scale_color_gradient2(midpoint= 0, high = "#CC6600", mid = "white", low = "#3399CC")+
  scale_fill_gradient2(midpoint= 0, high = "#CC6600", mid = "white", low = "#3399CC")+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )


corrPlot <- ggplot( data = Correlations, aes( y = Correlation, x = parameters ))+
  geom_jitter(shape = 21)+
  scale_y_continuous(expression("Kendall's "*tau*" (Phenotype, Optimum)"), limits = c(0,1))+
  scale_x_discrete("")+
  theme_bw()+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )


combinedPlot<- ggarrange(dataPlot, corrPlot, widths = c(3,1))

png("~/work/GEA/simulations/Plots/PhenEnvCorr.png", width = 10, height = 6, units="in", res = 300)
print(combinedPlot)
dev.off()

## Now add a couple of Manhattan plots for good measure...

gea_1 <- read.csv("~/work/GEA/simulations/G.0.2/88.checkMeOutChump.csv")
gea_2 <- read.csv("~/work/GEA/simulations/G.0.2/87.checkMeOutChump.csv")
gea_3 <- read.csv("~/work/GEA/simulations/G.0.2/88.checkMeOutChump.csv")
gea_4 <- read.csv("~/work/GEA/simulations/G.0.2/89.checkMeOutChump.csv")
gea_4$pve <- 2* gea_4$selCoeff^2 * gea_4$pbar_q_bar
`plot(gea_4$position, gea_4$pve)


geaPlot1 <- ggplot(data= gea_1, aes(x = position/1e6, y = -log10(pop_k_tau_p_value), col = as.factor(LA)))+
  geom_jitter( alpha = 0.3 , width = 0.0025)+
  #  geom_point( col = "#3399CC")+
  scale_y_continuous(expression("-log"[10]*"("*italic("p")*"-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  scale_color_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )

geaPlot2 <- ggplot(data= gea_2, aes(x = position/1e6, y = -log10(pop_k_tau_p_value), col = as.factor(LA)))+
  geom_jitter( alpha = 0.3 , width = 0.0025)+
  #  geom_point( col = "#3399CC")+
  scale_y_continuous(expression("-log"[10]*"("*italic("p")*"-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  scale_color_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )

geaPlot3 <- ggplot(data= gea_3, aes(x = position/1e6, y = -log10(pop_k_tau_p_value), col = as.factor(LA)))+
  geom_jitter( alpha = 0.3 , width = 0.0025)+
  #  geom_point( col = "#3399CC")+
  scale_y_continuous(expression("-log"[10]*"("*italic("p")*"-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  scale_color_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )
geaPlot4 <- ggplot(data= gea_4, aes(x = position/1e6, y = -log10(pop_k_tau_p_value), col = as.factor(LA)))+
  geom_jitter( alpha = 0.3 , width = 0.0025)+
  #  geom_point( col = "#3399CC")+
  scale_y_continuous(expression("-log"[10]*"("*italic("p")*"-value)"), limits = c(0,10))+
  scale_x_continuous("Position (Mbp)")+
  scale_color_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )
ggarrange(geaPlot1, geaPlot2, geaPlot3, nrow = 3, ncol =1 )




weiZ <- function(gea){
  gea$z <- qnorm(gea$pop_k_tau_p_value, lower.tail = F)
  gea$z[gea$z==-Inf] <- qnorm(0.999, lower.tail = F)

  gea_filt<-gea#[gea$maf > 0.05, ]

  gene_pos <- tapply( gea_filt$position, as.factor(gea_filt$gene), mean)
  
  LA <- tapply( gea_filt$LA, as.factor(gea_filt$gene), mean)

    SNPs <- tapply( gea_filt$LA, as.factor(gea_filt$gene), length)
  
  weiZ_num <- tapply( gea_filt$pbar_q_bar * gea_filt$z, as.factor(gea_filt$gene), sum)

  weiZ_den <- sqrt( tapply( gea_filt$pbar_q_bar**2, gea_filt$gene, sum))

  hits <- names( sort(weiZ_num/weiZ_den, decreasing = T)[1:5] )

#sort( names( summary( gea[gea$LA == 1,]$gene)) )

  weiDF <- data.frame( LA = LA, position = gene_pos,  weiZ =weiZ_num/ weiZ_den, nSNPs = SNPs )

  return(weiDF)
}

gea1_weiZ<- weiZ(gea_1)

gea1_weiZ_Plot <- ggplot(data = gea1_weiZ, aes( x = position/1e6, y = weiZ, col = as.factor(LA)))+
  geom_point(size = 2.5)+
  scale_y_continuous("WZA", limits = c(-2,10))+
  geom_hline( aes(yintercept = qnorm(0.05/50, lower.tail = F)), lty = 2, alpha = 0.5)+
  scale_x_continuous("Position (Mbp)")+
  scale_colour_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )

gea2_weiZ<- weiZ(gea_2)
gea2_weiZ_Plot <- ggplot(data = gea2_weiZ, aes( x = position/1e6, y = weiZ, col = as.factor(LA)))+
  geom_point(size = 2.5)+
  scale_y_continuous("WZA", limits = c(-2,10))+
  geom_hline( aes(yintercept = qnorm(0.05/50, lower.tail = F)), lty = 2, alpha = 0.5)+
  scale_x_continuous("Position (Mbp)")+
  scale_colour_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )

gea3_weiZ<- weiZ(gea_3)
gea3_weiZ_Plot <- ggplot(data = gea3_weiZ, aes( x = position/1e6, y = weiZ, col = as.factor(LA)))+
  geom_point(size = 2.5)+
  scale_y_continuous("WZA", limits = c(-2,10))+
  geom_hline( aes(yintercept = qnorm(0.05/50, lower.tail = F)), lty = 2, alpha = 0.5)+
  scale_x_continuous("Position (Mbp)")+
  scale_colour_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )


gea4_weiZ<- weiZ(gea_4)
gea4_weiZ_Plot <- ggplot(data = gea4_weiZ, aes( x = position/1e6, y = weiZ, col = as.factor(LA)))+
  geom_point(size = 2.5)+
  scale_y_continuous("WZA", limits = c(-2,10))+
  geom_hline( aes(yintercept = qnorm(0.05/50, lower.tail = F)), lty = 2, alpha = 0.5)+
  scale_x_continuous("Position (Mbp)")+
  scale_colour_manual(values = c("#3399CC","#CC6600"))+
  theme_bw()+
  guides(colour=F, fill=F)+
  theme(
    axis.title.x =  element_text(size = 14, angle= 0, colour = "black"),
    axis.title.y =  element_text(size = 14, colour = "black"),
    axis.text.x =  element_text(size = 12, angle= 0, colour = "black"),
    axis.text.y =  element_text(size = 12, angle= 0, colour = "black")
  )



png("~/work/GEA/simulations/Plots/ManhattanPlot.png", width = 8, height = 10, units="in", res = 300)
ggarrange(geaPlot1, gea1_weiZ_Plot, geaPlot2, gea2_weiZ_Plot, geaPlot3, gea3_weiZ_Plot,geaPlot4, gea4_weiZ_Plot, nrow = 4, ncol =2 )
dev.off()

