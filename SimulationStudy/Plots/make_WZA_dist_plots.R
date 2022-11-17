rm(list = ls())

weak_bc <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/BC_Map_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")
weak_gradient <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/cline_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")
weak_trunc <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/trunc_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")

weak_bc$map <- "BC Map"
weak_gradient$map <- "Gradient Map"
weak_trunc$map <- "Truncated Map"

weak <- rbind( weak_bc, weak_gradient, weak_trunc)
weak$selection <- "Weak Selection"

strong_bc <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/BC_Map_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")
strong_gradient <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/cline_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")
strong_trunc <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/trunc_output_n40_i20/1_WZA_geno_k_tau_p_value.csv")

strong_bc$map <- "BC Map"
strong_gradient$map <- "Gradient Map"
strong_trunc$map <- "Truncated Map"

strong <- rbind( strong_bc, strong_gradient, strong_trunc)
strong$selection <- "Strong Selection"

all <- rbind(strong, weak)

library(ggplot2)

WZA_dist <- ggplot(data = all, aes(x = Z, fill = LA>0.005))+
  geom_density(bw=1, alpha = 0.5)+
  facet_grid(selection~map)+
  xlab(expression(italic(Z[W])))+
  ylab("Density")+
  scale_fill_manual("Locally Adaptive",
                    values = c("grey","#56B4E9") )+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("~/work/GEA/simulations/SimulationStudy/Plots/WZA_distribution.pdf",
       width = 10, height = 5)
