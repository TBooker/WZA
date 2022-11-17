rm(list = ls())

## Make the distribution of effect size plots

bc_strong <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/plotting_files/BC_Map_WZA.out")
bc_strong$map <- "BC Map"
cline_strong <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/plotting_files/cline_WZA.out") 
cline_strong$map <- "Gradient Map"
trunc_strong  <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/plotting_files/trunc_WZA.out")
trunc_strong$map <- "Truncated Map"
strong <- rbind( bc_strong, cline_strong, trunc_strong)
strong$effect <- "Strong Selection"

bc_weak  <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/plotting_files/BC_Map_WZA.out")
bc_weak$map <- "BC Map"
cline_weak  <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/plotting_files/cline_WZA.out")
cline_weak$map <- "Gradient Map"
trunc_weak  <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/plotting_files/trunc_WZA.out")
trunc_weak$map <- "Truncated Map"
weak <- rbind( bc_weak, cline_weak, trunc_weak)
weak$effect <- "Weak Selection"

all <- rbind(weak, strong)
all_screend <- all[all$LA!=0,]
all_screend <- all_screend[all_screend$LA!=-99,]

library(ggplot2)

effect_size_graph <- ggplot(data = all_screend, 
       aes(x = LA,
           fill = effect))+
  geom_histogram(bins = 40, col="black")+
  xlab("Measure of Locally Adaptive Effect")+
  ylab("Count")+
  geom_vline(aes(xintercept = 0.005), lty = 2, col = "red", lwd = 1)+
  scale_fill_manual(values = c("#1B9E77", "#D95F02" ))+
  facet_grid(effect~map)+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13),
        legend.title = element_blank())


ggsave(filename = "~/work/GEA/simulations/SimulationStudy/Plots/effect_size_graph.pdf", 
       effect_size_graph,
       width = 12,
       height = 4)
