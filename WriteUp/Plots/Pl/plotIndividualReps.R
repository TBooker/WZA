rm(list=ls())



library(ggplot2)
library(ggpubr)


all <- read.csv("/media/booker/HOWDY/GEA/G.1.1/cline_d40_n50_WZA.csv")

for (Vs in c( "100")){
  for (Ua in c("0.001", "0.0001")){
    for (sig_a in c("0.3", "1")){
      pdf(paste("~/work/GEA/inspectionPlots/cline_map_",sig_a,"_",Ua,"_",Vs,".pdf", sep = ''), height = 10, width = 6)
      for (rep in 21:40){
      
        rep1 <- read.csv( paste("/media/booker/HOWDY/GEA/G.1.1/cline/",rep,"_",sig_a,"_",Ua,"_",Vs,".csv" , sep = ''))
        
        A<- ggplot(data = rep1, aes( x= position/1e6, y = -log10(geno_k_tau_p_value), col = LA > 0.01))+
          geom_point(alpha = 0.67)+
          ggtitle( paste("Replicate =",rep))+
          scale_color_manual( values = c("#009999","#FF7400"))+
          scale_x_continuous("Position (Mbp)", limits = c(0,1))+
          scale_y_continuous("-log10(p-value)")+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5))
        
        temp_slice = all[(all$U_a == as.numeric(Ua))&(all$sig_a == as.numeric(sig_a))&(all$Vs == as.numeric(Vs))&(all$rep == rep),]
        
        B<- ggplot(data = temp_slice, aes( x= position/1e6, y = Z, col = LA>0.01 ))+
          geom_point(size = 2.5)+
          scale_color_manual( values = c("#009999","#FF7400"))+
          scale_x_continuous("Position (Mbp)", limits = c(0,1))+
          scale_y_continuous(limits = c(0,30))+
          theme_bw()
        
        C<- ggplot(data = temp_slice, aes( x= position/1e6, y = -log10(top_candidate_p), col = LA > 0.01))+
          geom_point(size = 2.5)+
          scale_color_manual( values = c("#009999","#FF7400"))+
          scale_x_continuous("Position (Mbp)", limits = c(0,1))+
          scale_y_continuous("-log10(top-candidate p-value)", limits = c(0,15))+
          theme_bw()
        
        print( ggarrange( A, B, C, ncol = 1, nrow = 3, align = "v" ) )
      }
    
      dev.off()
    }
  }
}
