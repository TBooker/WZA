rm( list = ls() )

## Here's a script for  
cline_directional <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/cline_sampled.WZA.csv")
BC_directional <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/BC_Map_sampled.WZA.csv")
trunc_directional <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/trunc_sampled.WZA.csv")

cline_directional$selection <- "Directional selection"
BC_directional$selection <- "Directional selection"
trunc_directional$selection <- "Directional selection"

cline_directional$map <- "Gradient map"
BC_directional$map <- "BC map"
trunc_directional$map <- "Truncated map"

cline_stabilising <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_Vs192/cline_sampled.WZA.csv")
BC_stabilising <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_Vs192/BC_Map_sampled.WZA.csv")
trunc_stabilising <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_Vs192/trunc_sampled.WZA.csv")

cline_stabilising$selection <- "Stabilising selection"
BC_stabilising$selection <- "Stabilising selection"
trunc_stabilising$selection <- "Stabilising selection"

cline_stabilising$map <- "Gradient map"
BC_stabilising$map <- "BC map"
trunc_stabilising$map <- "Truncated map"


library(ggplot2)

Z_data <- rbind( cline_stabilising,
                 BC_stabilising,
                 trunc_stabilising,
                 cline_directional,
                 BC_directional,
                 trunc_directional)

selectionWZA_plot<-ggplot(data = Z_data[Z_data$LA>-1,], aes( x = empR_Z, fill = LA>0.005))+
  geom_density( alpha = 0.75 )+
  facet_grid(selection ~ map)+
  scale_x_continuous(expression(italic(Z[W])))+
  scale_y_continuous("Density", limits = c(0,0.5))+
#  stat_summary(aes(x = 0.1, y = empR_Z, xintercept = stat(y), col = LA>0.005), 
#               fun = mean, geom = "vline") +
  scale_fill_manual("Locally adaptive gene", 
                    values = c("black","lightblue"))+
  theme_bw()+
  theme(
    strip.text = element_text( face = "italic" ),
    strip.background = element_blank()
  )

pdf("~/UBC/GEA/WZA/WriteUp/Plots/selectionWZA_plot.pdf", width= 10, height = 6)
print(selectionWZA_plot)
dev.off()


## How many adapted genes are there in any simulation?


mean(tapply( BC_directional$LA, BC_directional$rep, function(x) sum(x>0.005)))
mean(tapply( cline_directional$LA, cline_directional$rep, function(x) sum(x>0.005)))
mean(tapply( trunc_directional$LA, trunc_directional$rep, function(x) sum(x>0.005)))



mean(tapply( BC_stabilising$LA, BC_stabilising$rep, function(x) sum(x>0.005)))
mean(tapply( cline_stabilising$LA, cline_stabilising$rep, function(x) sum(x>0.005)))
mean(tapply( trunc_stabilising$LA, trunc_stabilising$rep, function(x) sum(x>0.005)))
