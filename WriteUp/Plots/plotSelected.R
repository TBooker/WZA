rm(list= ls())
library(ggplot2)
library(ggpubr)


## Let's start by making plots for the unsampled data
set_0.3_0.001_100_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.001_100_all.csv")
set_0.3_0.0001_100_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.0001_100_all.csv")
set_0.3_0.001_200_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.001_200_all.csv")
set_0.3_0.0001_200_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.0001_200_all.csv")
set_1_0.001_100_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.001_100_all.csv")
set_1_0.0001_100_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.0001_100_all.csv")
set_1_0.001_200_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.001_200_all.csv")
set_1_0.0001_200_all<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.0001_200_all.csv")

allDemes<-read.csv("/media/booker/HOWDY/GEA/G.1.0/")

#allDemes<-rbind(set_0.3_0.0001_100_all,set_0.3_0.0001_200_all,set_0.3_0.001_100_all,set_0.3_0.001_200_all,set_1_0.0001_100_all,set_1_0.0001_200_all,set_1_0.001_100_all,set_1_0.001_200_all) 
allDemes$LA2<- allDemes$LA>0.005
  
allDemes$sig_a <- factor(as.factor(allDemes$sig_a),
                    levels = c("0.3",
                               "1"), 
                    labels = c(expression(italic(sigma[a]^"2")*" = 0.3"),
                               expression(italic(sigma[a]^"2")*" = 1.0" )
                    ))

allDemes$U_a <- factor(as.factor(allDemes$U_a),
                  levels = c("0.001",
                             "1e-04"), 
                  labels = c(expression(italic(U[a])*" = 0.001"),
                             expression(italic(U[a])*" = 0.0001" )
                  ))

allDemes$Vs <- factor(as.factor(allDemes$Vs),
                 levels = c("100",
                            "200"), 
                 labels = c(expression(V[s]*" = 100"),
                            expression(V[s]*" = 200" )
                 ))

## PVE Vs. Zw

pdf("~/work/GEA/simulations/Plots/PVE_v_Z.pdf", width = 8, height = 8)

ggplot(data = allDemes[allDemes$LA > 0.05,], aes(x= LA, y = Z))+
  geom_point(alpha = 0.5)+
  facet_grid(  sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("PVE")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  guides(color = F)+
  scale_y_continuous(expression(italic(Z[w])))+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )

dev.off()



## Weighted Z test

pdf("~/work/GEA/simulations/Plots/weightedZ_all.pdf", width = 8, height = 8)
ggplot(data = allDemes[allDemes$LA2 == F,], aes(x= SNPs, y = Z, col = LA2))+
  geom_point(alpha = 0.5)+
  geom_point(data = allDemes[allDemes$LA2 == T,], aes(x= SNPs, y = Z, col = LA2), alpha = 0.8)+
  facet_grid( sig_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("Number of SNPs")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  scale_y_continuous(expression(italic(Z[w])))+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )
dev.off()


pdf("~/work/GEA/simulations/Plots/topCan_all.pdf", width = 8, height = 8)
ggplot(data = allDemes[allDemes$LA2 == F,], aes(x= SNPs, y = hits, col = LA2))+
  geom_point(alpha = 0.5)+
  geom_point(data = allDemes[allDemes$LA2 == T,], aes(x= SNPs, y = hits, col = LA2), alpha = 0.8)+
  facet_grid( sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("Number of SNPs")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  scale_y_continuous("Number of Hits")+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )
dev.off()


## Let's start by making plots for the sampled data
set_0.3_0.001_100_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.001_100_d40.csv")
set_0.3_0.0001_100_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.0001_100_d40.csv")
set_0.3_0.001_200_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.001_200_d40.csv")
set_0.3_0.0001_200_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/0.3_0.0001_200_d40.csv")
set_1_0.001_100_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.001_100_d40.csv")
set_1_0.0001_100_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.0001_100_d40.csv")
set_1_0.001_200_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.001_200_d40.csv")
set_1_0.0001_200_d40<-read.csv("/media/booker/HOWDY/GEA/G.0.8_selected2_analysis/1_0.0001_200_d40.csv")

sampledDemes<-rbind(set_0.3_0.0001_100_d40,set_0.3_0.0001_200_d40,set_0.3_0.001_100_d40,set_0.3_0.001_200_d40,set_1_0.0001_100_d40,set_1_0.0001_200_d40,set_1_0.001_100_d40,set_1_0.001_200_d40) 
sampledDemes<-read.csv("/media/booker/HOWDY/GEA/G.1.0/BC_Map.WZA.d40_n50.csv")
sampledDemes$LA2<- sampledDemes$LA>0.005

sampledDemes$sig_a <- factor(as.factor(sampledDemes$sig_a),
                         levels = c("0.3",
                                    "1"), 
                         labels = c(expression(italic(sigma[a]^"2")*" = 0.3"),
                                    expression(italic(sigma[a]^"2")*" = 1.0" )
                         ))

sampledDemes$U_a <- factor(as.factor(sampledDemes$U_a),
                       levels = c("0.001",
                                  "1e-04"), 
                       labels = c(expression(italic(U[a])*" = 0.001"),
                                  expression(italic(U[a])*" = 0.0001" )
                       ))

sampledDemes$Vs <- factor(as.factor(sampledDemes$Vs),
                      levels = c("100",
                                 "200"), 
                      labels = c(expression(V[s]*" = 100"),
                                 expression(V[s]*" = 200" )
                      ))



ggplot(data = sampledDemes[sampledDemes$LA > 0.005,], aes(x= LA, y = Z))+
  geom_point(alpha = 0.5)+
  facet_grid(  sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("PVE")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  guides(color = F)+
  scale_y_continuous(expression(italic(Z[w])))+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )

## Weighted Z test

#pdf("~/work/GEA/simulations/Plots/weightedZ_sample.pdf", width = 8, height = 8)
ggplot(data = sampledDemes[sampledDemes$LA2 == F,], aes(x= SNPs, y = Z, col = LA2))+
  geom_point(alpha = 0.5)+
  geom_point(data = sampledDemes[sampledDemes$LA2 == T,], aes(x= SNPs, y = Z, col = LA2), alpha = 0.8)+
  facet_grid(   sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("Number of SNPs")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  scale_y_continuous(expression(italic(Z[w])))+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )
#dev.off()

sigLine_SNPs <- 10:75
sigLine_hits <- qbinom(.9999,sigLine_SNPs,0.01)
lineDF<- data.frame( SNPs = sigLine_SNPs, hits = sigLine_hits)

pdf("~/work/GEA/simulations/Plots/topCan_sample.pdf", width = 8, height = 8)

ggplot(data = sampledDemes[sampledDemes$LA2 == F,], aes(x= SNPs, y = -log10(top_candidate_p), col = LA2))+
  geom_point(alpha = 0.5)+
  geom_point(data = sampledDemes[sampledDemes$LA2 == T,], aes(x= SNPs, y = hits, col = LA2), alpha = 0.8)+
  geom_line(data = lineDF, aes( x = SNPs, y = hits), col = "black")+
  facet_grid(   sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("Number of SNPs")+
  scale_color_brewer(">0.005 PVE", palette = "Set2")+
  scale_y_continuous("Number of Hits")+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )
dev.off()


allDemes<-read.csv("/media/booker/HOWDY/GEA/G.1.0/XXX")
hist(allDemes$Z)

allDemes<-allDemes[(allDemes$sig_a == 0.3)&(allDemes$Vs == 200)&(allDemes$U_a == 0.0001),]
allDemes$sigZ <-allDemes$Z>25
allDemes$sigTC <-allDemes$top_candidate_p<0.0001 
allDemes$sigN <-allDemes$bigHits>0 
model <- glm(sigZ~LA, family = binomial(link = "logit"), data = allDemes[allDemes$LA > 0.05,])
model2 <- glm(sigTC~LA, family = binomial(link = "logit"), data = allDemes[allDemes$LA > 0,])
model3 <- glm(sigN~LA, family = binomial(link = "logit"), data = allDemes[allDemes$LA > 0,])


plot(allDemes$sigTC~allDemes$LA)
xweight <- seq(0,1,0.001)
yweight <- predict(model, list(LA = xweight), type = "response")
yweight2 <- predict(model2, list(LA = xweight), type = "response")
yweight3 <- predict(model3, list(LA = xweight), type = "response")
lines(xweight, yweight)
lines(xweight, yweight2, col = "red")
lines(xweight, yweight3, col = "blue")

ggplot(data = allDemes[allDemes$LA > 0.05,], aes(x= LA, y = Z))+
  geom_point(alpha = 0.5)+
  facet_grid(  sig_a + U_a ~ Vs, labeller = label_parsed)+
  scale_x_continuous("PVE")+
  scale_color_brewer(">0.05 PVE", palette = "Set2")+
  guides(color = F)+
  scale_y_continuous(expression(italic(Z[w])))+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )

