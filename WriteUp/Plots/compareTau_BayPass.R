rm(list = ls())

bp<-read.table("/media/booker/HOWDY/GEA/G.1.1/BC_map/6_1_0.001_100_d40_n50_summary_betai_reg.out", header = T)
csv<-read.csv("/media/booker/HOWDY/GEA/G.1.1/BC_map/6_1_0.001_100_d40_n50.csv",)
csv$bf<- bp$BF.dB.

library(ggplot2)
library(ggpubr)

ggplot( data = csv, aes(x = position, y = -log10(pop_k_tau_p_value), col = as.factor(LA > 0.001) ))+
  geom_point()
ggplot( data = csv, aes(x = position, y = bf, col = as.factor(LA > 0.001) ))+
  geom_point()
ggplot( data = csv, aes(x = -log10(pop_k_tau_p_value), y = bf, col = as.factor(LA > 0.001) ))+
  geom_point(alpha= 0.75)

cor(-log10(csv$pop_k_tau_p_value), csv$bf )



Z<-read.csv("/media/booker/HOWDY/GEA/G.1.1/BC_map_d40_n50_WZA_bayPass.csv")

Z$sig_a <- factor(Z$sig_a ,levels = c(0.3,1),
                               labels = c(expression(sigma[alpha]^2*" = 0.3" ),expression(sigma[alpha]^2*" = 1.0" )))

Z$U_a <- factor(Z$U_a , levels = c(0.001,1e-04),
                            labels = c(expression(U[alpha]*" = 0.001" ),expression(U[alpha]*" = 0.0001" )))


Z$Vs <- factor(Z$Vs , levels = c(100,200),
                            labels = c(expression(V[s]*" = 100" ),expression(V[s]*" = 200" )))

Z$POS = Z$position + (Z$rep*1000000)
Z_ken <- ggplot( data = Z, aes(x = POS/1e6, y = Z_kendall, col = as.factor(LA > 0.05) ))+
  geom_point()+
  facet_wrap(U_a+sig_a+Vs~., nrow = 1, labeller = label_parsed)+
  scale_y_continuous(expression("Weighted Z"["Kendall's"]))+
    scale_color_manual(values = c("black","red"))+
  theme_bw()

TC <- ggplot( data = Z, aes(x = POS/1e6, y = -log10(top_candidate_p), col = as.factor(LA > 0.05) ))+
  geom_point()+
  facet_wrap(U_a+sig_a+Vs~., nrow = 1, labeller = label_parsed)+
  scale_y_continuous(expression("Top Candidate Test"["Kendall's"]))+
  scale_color_manual(values = c("black","red"))+
  theme_bw()

Z_bp <- ggplot( data = Z, aes(x = POS/1e6, y = Z_bayPass, col = as.factor(LA > 0.05) ))+
  geom_point()+
  facet_wrap(U_a+sig_a+Vs~., nrow = 1, labeller = label_parsed)+
  scale_y_continuous(expression("Weighted Z"["BayPass"]))+
  scale_color_manual(values = c("black","red"))+
    theme_bw()

png("/media/booker/HOWDY/GEA/Plots/BC_Map_statComparison.png", height = 8, width = 14, res = 200, units = "in")
ggarrange(Z_ken, TC, Z_bp, nrow = 3, ncol =1, common.legend = T)
dev.off()

