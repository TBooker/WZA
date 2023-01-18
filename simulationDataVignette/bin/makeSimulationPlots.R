rm(list = ls())

## Make some plots for the GitHub simulation vignette...
# single SNP GEA...

gea <- read.csv("../output/simulationSingleSNP_GEA.txt",
         sep = "\t")

library(ggplot2)

gea<-gea[gea$MAF>0.05,]

sig_thresh <- 1-(1-0.05)^(1/nrow(gea))

gea_point <-ggplot(data = gea, 
       aes(x = position/1e6, 
           y = -log10(K_tau_p),
           fill = LA))+
  geom_point(shape = 21,
             size = 2)+
  geom_hline(yintercept=-log10(sig_thresh),
             lty = 2)+
  scale_fill_manual("Locally Adaptive Gene",
                    values = c("white","black"))+
  ylab(expression("-log"[10]*"("*italic(p)*"-value)"))+
  xlab("Position in Genome (Mbp)")+
  theme_bw()

ggsave(filename = "../plots/GEA_points.png",
       plot = gea_point,
       height = 4,
       width = 10)


gea_hist <-ggplot(data = gea, 
                   aes(x = K_tau_p))+
  geom_histogram(col = "black",
                 fill = "white",
                 bins = 50)+
  xlab(expression("-log"[10]*"("*italic(p)*"-value)"))+
  ylab("Count")+
  theme_bw()
              
ggsave(filename = "../plots/GEA_hist.png",
       plot = gea_hist,
       height = 5,
       width = 8)





# Running the WZA

wza <- read.csv("../output/simulation_wza.csv")

sig_thresh_wza <- 1-(1-0.05)^(1/nrow(wza))

wza_point <-ggplot(data = wza, 
                   aes(x = position/1e6, 
                       y = -log10(Z_pVal),
                       fill = LA==1))+
  geom_point(shape = 21,
             size = 2)+
  geom_hline(yintercept=-log10(sig_thresh_wza),
             lty = 2)+
  scale_fill_manual("Locally Adaptive Gene",
                    values = c("white","black"))+
  ylab(expression("-log"[10]*"("*italic(p)*"-value)"))+
  xlab("Position in Genome (Mbp)")+
  theme_bw()
  

ggsave(filename = "../plots/WZA_points.png",
       plot = wza_point,
       height = 4,
       width = 10)
  
  
wza_hist <-ggplot(data = wza, 
                    aes(x = Z_pVal))+
    geom_histogram(col = "black",
                   fill = "white",
                   bins = 30)+
    xlab(expression("-log"[10]*"("*italic(p)*"-value)"))+
    ylab("Count")+
    theme_bw()

ggsave(filename = "../plots/WZA_hist.png",
       plot = wza_hist,
       height = 5,
       width = 8)
