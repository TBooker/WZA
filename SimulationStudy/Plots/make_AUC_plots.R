rm(list = ls())
library(ggplot2)
library(reshape2)
##################### BC Map - uniform

s0.0136_bc_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i5")
s0.0136_bc_n40_i5$demes <- "40 Demes"
s0.0136_bc_n40_i5$inds <- "5 Individuals"

s0.0136_bc_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i10")
s0.0136_bc_n40_i10$demes <- "40 Demes"
s0.0136_bc_n40_i10$inds <- "10 Individuals"

s0.0136_bc_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20")
s0.0136_bc_n40_i20$demes <- "40 Demes"
s0.0136_bc_n40_i20$inds <- "20 Individuals"

s0.0136_bc_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i5")
s0.0136_bc_n20_i5$demes <- "20 Demes"
s0.0136_bc_n20_i5$inds <- "5 Individuals"

s0.0136_bc_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i10")
s0.0136_bc_n20_i10$demes <- "20 Demes"
s0.0136_bc_n20_i10$inds <- "10 Individuals"

s0.0136_bc_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20")
s0.0136_bc_n20_i20$demes <- "20 Demes"
s0.0136_bc_n20_i20$inds <- "20 Individuals"

s0.0136_bc_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i5")
s0.0136_bc_n10_i5$demes <- "10 Demes"
s0.0136_bc_n10_i5$inds <- "5 Individuals"

s0.0136_bc_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i10")
s0.0136_bc_n10_i10$demes <- "10 Demes"
s0.0136_bc_n10_i10$inds <- "10 Individuals"

s0.0136_bc_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20")
s0.0136_bc_n10_i20$demes <- "10 Demes"
s0.0136_bc_n10_i20$inds <- "20 Individuals"

s0.0136_bc <-rbind( s0.0136_bc_n40_i5,
       s0.0136_bc_n40_i10,
       s0.0136_bc_n40_i20,
       s0.0136_bc_n20_i5,
       s0.0136_bc_n20_i10,
       s0.0136_bc_n20_i20,
       s0.0136_bc_n10_i5,
       s0.0136_bc_n10_i10,
       s0.0136_bc_n10_i20)



s0.003_bc_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i5")
s0.003_bc_n40_i5$demes <- "40 Demes"
s0.003_bc_n40_i5$inds <- "5 Individuals"

s0.003_bc_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i10")
s0.003_bc_n40_i10$demes <- "40 Demes"
s0.003_bc_n40_i10$inds <- "10 Individuals"

s0.003_bc_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n40_i20")
s0.003_bc_n40_i20$demes <- "40 Demes"
s0.003_bc_n40_i20$inds <- "20 Individuals"

s0.003_bc_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i5")
s0.003_bc_n20_i5$demes <- "20 Demes"
s0.003_bc_n20_i5$inds <- "5 Individuals"

s0.003_bc_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i10")
s0.003_bc_n20_i10$demes <- "20 Demes"
s0.003_bc_n20_i10$inds <- "10 Individuals"

s0.003_bc_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n20_i20")
s0.003_bc_n20_i20$demes <- "20 Demes"
s0.003_bc_n20_i20$inds <- "20 Individuals"

s0.003_bc_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i5")
s0.003_bc_n10_i5$demes <- "10 Demes"
s0.003_bc_n10_i5$inds <- "5 Individuals"

s0.003_bc_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i10")
s0.003_bc_n10_i10$demes <- "10 Demes"
s0.003_bc_n10_i10$inds <- "10 Individuals"

s0.003_bc_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_n10_i20")
s0.003_bc_n10_i20$demes <- "10 Demes"
s0.003_bc_n10_i20$inds <- "20 Individuals"

s0.003_bc <-rbind( s0.003_bc_n40_i5,
                    s0.003_bc_n40_i10,
                    s0.003_bc_n40_i20,
                    s0.003_bc_n20_i5,
                    s0.003_bc_n20_i10,
                    s0.003_bc_n20_i20,
                    s0.003_bc_n10_i5,
                    s0.003_bc_n10_i10,
                    s0.003_bc_n10_i20)



##################### cline - uniform

s0.0136_cline_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i5")
s0.0136_cline_n40_i5$demes <- "40 Demes"
s0.0136_cline_n40_i5$inds <- "5 Individuals"

s0.0136_cline_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i10")
s0.0136_cline_n40_i10$demes <- "40 Demes"
s0.0136_cline_n40_i10$inds <- "10 Individuals"

s0.0136_cline_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20")
s0.0136_cline_n40_i20$demes <- "40 Demes"
s0.0136_cline_n40_i20$inds <- "20 Individuals"

s0.0136_cline_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i5")
s0.0136_cline_n20_i5$demes <- "20 Demes"
s0.0136_cline_n20_i5$inds <- "5 Individuals"

s0.0136_cline_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i10")
s0.0136_cline_n20_i10$demes <- "20 Demes"
s0.0136_cline_n20_i10$inds <- "10 Individuals"

s0.0136_cline_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20")
s0.0136_cline_n20_i20$demes <- "20 Demes"
s0.0136_cline_n20_i20$inds <- "20 Individuals"

s0.0136_cline_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i5")
s0.0136_cline_n10_i5$demes <- "10 Demes"
s0.0136_cline_n10_i5$inds <- "5 Individuals"

s0.0136_cline_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i10")
s0.0136_cline_n10_i10$demes <- "10 Demes"
s0.0136_cline_n10_i10$inds <- "10 Individuals"

s0.0136_cline_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20")
s0.0136_cline_n10_i20$demes <- "10 Demes"
s0.0136_cline_n10_i20$inds <- "20 Individuals"

s0.0136_cline <-rbind( s0.0136_cline_n40_i5,
                    s0.0136_cline_n40_i10,
                    s0.0136_cline_n40_i20,
                    s0.0136_cline_n20_i5,
                    s0.0136_cline_n20_i10,
                    s0.0136_cline_n20_i20,
                    s0.0136_cline_n10_i5,
                    s0.0136_cline_n10_i10,
                    s0.0136_cline_n10_i20)



s0.003_cline_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i5")
s0.003_cline_n40_i5$demes <- "40 Demes"
s0.003_cline_n40_i5$inds <- "5 Individuals"

s0.003_cline_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i10")
s0.003_cline_n40_i10$demes <- "40 Demes"
s0.003_cline_n40_i10$inds <- "10 Individuals"

s0.003_cline_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n40_i20")
s0.003_cline_n40_i20$demes <- "40 Demes"
s0.003_cline_n40_i20$inds <- "20 Individuals"

s0.003_cline_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i5")
s0.003_cline_n20_i5$demes <- "20 Demes"
s0.003_cline_n20_i5$inds <- "5 Individuals"

s0.003_cline_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i10")
s0.003_cline_n20_i10$demes <- "20 Demes"
s0.003_cline_n20_i10$inds <- "10 Individuals"

s0.003_cline_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n20_i20")
s0.003_cline_n20_i20$demes <- "20 Demes"
s0.003_cline_n20_i20$inds <- "20 Individuals"

s0.003_cline_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i5")
s0.003_cline_n10_i5$demes <- "10 Demes"
s0.003_cline_n10_i5$inds <- "5 Individuals"

s0.003_cline_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i10")
s0.003_cline_n10_i10$demes <- "10 Demes"
s0.003_cline_n10_i10$inds <- "10 Individuals"

s0.003_cline_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_n10_i20")
s0.003_cline_n10_i20$demes <- "10 Demes"
s0.003_cline_n10_i20$inds <- "20 Individuals"

s0.003_cline <-rbind( s0.003_cline_n40_i5,
                   s0.003_cline_n40_i10,
                   s0.003_cline_n40_i20,
                   s0.003_cline_n20_i5,
                   s0.003_cline_n20_i10,
                   s0.003_cline_n20_i20,
                   s0.003_cline_n10_i5,
                   s0.003_cline_n10_i10,
                   s0.003_cline_n10_i20)




##################### trunc - uniform

s0.0136_trunc_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i5")
s0.0136_trunc_n40_i5$demes <- "40 Demes"
s0.0136_trunc_n40_i5$inds <- "5 Individuals"

s0.0136_trunc_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i10")
s0.0136_trunc_n40_i10$demes <- "40 Demes"
s0.0136_trunc_n40_i10$inds <- "10 Individuals"

s0.0136_trunc_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20")
s0.0136_trunc_n40_i20$demes <- "40 Demes"
s0.0136_trunc_n40_i20$inds <- "20 Individuals"

s0.0136_trunc_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i5")
s0.0136_trunc_n20_i5$demes <- "20 Demes"
s0.0136_trunc_n20_i5$inds <- "5 Individuals"

s0.0136_trunc_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i10")
s0.0136_trunc_n20_i10$demes <- "20 Demes"
s0.0136_trunc_n20_i10$inds <- "10 Individuals"

s0.0136_trunc_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20")
s0.0136_trunc_n20_i20$demes <- "20 Demes"
s0.0136_trunc_n20_i20$inds <- "20 Individuals"

s0.0136_trunc_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i5")
s0.0136_trunc_n10_i5$demes <- "10 Demes"
s0.0136_trunc_n10_i5$inds <- "5 Individuals"

s0.0136_trunc_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i10")
s0.0136_trunc_n10_i10$demes <- "10 Demes"
s0.0136_trunc_n10_i10$inds <- "10 Individuals"

s0.0136_trunc_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20")
s0.0136_trunc_n10_i20$demes <- "10 Demes"
s0.0136_trunc_n10_i20$inds <- "20 Individuals"

s0.0136_trunc <-rbind( s0.0136_trunc_n40_i5,
                        s0.0136_trunc_n40_i10,
                        s0.0136_trunc_n40_i20,
                        s0.0136_trunc_n20_i5,
                        s0.0136_trunc_n20_i10,
                        s0.0136_trunc_n20_i20,
                        s0.0136_trunc_n10_i5,
                        s0.0136_trunc_n10_i10,
                        s0.0136_trunc_n10_i20)


s0.003_trunc_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i5")
s0.003_trunc_n40_i5$demes <- "40 Demes"
s0.003_trunc_n40_i5$inds <- "5 Individuals"

s0.003_trunc_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i10")
s0.003_trunc_n40_i10$demes <- "40 Demes"
s0.003_trunc_n40_i10$inds <- "10 Individuals"

s0.003_trunc_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n40_i20")
s0.003_trunc_n40_i20$demes <- "40 Demes"
s0.003_trunc_n40_i20$inds <- "20 Individuals"

s0.003_trunc_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i5")
s0.003_trunc_n20_i5$demes <- "20 Demes"
s0.003_trunc_n20_i5$inds <- "5 Individuals"

s0.003_trunc_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i10")
s0.003_trunc_n20_i10$demes <- "20 Demes"
s0.003_trunc_n20_i10$inds <- "10 Individuals"

s0.003_trunc_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n20_i20")
s0.003_trunc_n20_i20$demes <- "20 Demes"
s0.003_trunc_n20_i20$inds <- "20 Individuals"

s0.003_trunc_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i5")
s0.003_trunc_n10_i5$demes <- "10 Demes"
s0.003_trunc_n10_i5$inds <- "5 Individuals"

s0.003_trunc_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i10")
s0.003_trunc_n10_i10$demes <- "10 Demes"
s0.003_trunc_n10_i10$inds <- "10 Individuals"

s0.003_trunc_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_n10_i20")
s0.003_trunc_n10_i20$demes <- "10 Demes"
s0.003_trunc_n10_i20$inds <- "20 Individuals"

s0.003_trunc <-rbind( s0.003_trunc_n40_i5,
                       s0.003_trunc_n40_i10,
                       s0.003_trunc_n40_i20,
                       s0.003_trunc_n20_i5,
                       s0.003_trunc_n20_i10,
                       s0.003_trunc_n20_i20,
                       s0.003_trunc_n10_i5,
                       s0.003_trunc_n10_i10,
                       s0.003_trunc_n10_i20)


#######################################################
#######################################################

##################### BC Map - mutMap

s0.0136_bc_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i5")
s0.0136_bc_mutMap_n40_i5$demes <- "40 Demes"
s0.0136_bc_mutMap_n40_i5$inds <- "5 Individuals"

s0.0136_bc_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i10")
s0.0136_bc_mutMap_n40_i10$demes <- "40 Demes"
s0.0136_bc_mutMap_n40_i10$inds <- "10 Individuals"

s0.0136_bc_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i20")
s0.0136_bc_mutMap_n40_i20$demes <- "40 Demes"
s0.0136_bc_mutMap_n40_i20$inds <- "20 Individuals"

s0.0136_bc_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i5")
s0.0136_bc_mutMap_n20_i5$demes <- "20 Demes"
s0.0136_bc_mutMap_n20_i5$inds <- "5 Individuals"

s0.0136_bc_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i10")
s0.0136_bc_mutMap_n20_i10$demes <- "20 Demes"
s0.0136_bc_mutMap_n20_i10$inds <- "10 Individuals"

s0.0136_bc_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i20")
s0.0136_bc_mutMap_n20_i20$demes <- "20 Demes"
s0.0136_bc_mutMap_n20_i20$inds <- "20 Individuals"

s0.0136_bc_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i5")
s0.0136_bc_mutMap_n10_i5$demes <- "10 Demes"
s0.0136_bc_mutMap_n10_i5$inds <- "5 Individuals"

s0.0136_bc_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i10")
s0.0136_bc_mutMap_n10_i10$demes <- "10 Demes"
s0.0136_bc_mutMap_n10_i10$inds <- "10 Individuals"

s0.0136_bc_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i20")
s0.0136_bc_mutMap_n10_i20$demes <- "10 Demes"
s0.0136_bc_mutMap_n10_i20$inds <- "20 Individuals"

s0.0136_bc_mutMap <-rbind( s0.0136_bc_mutMap_n40_i5,
                    s0.0136_bc_mutMap_n40_i10,
                    s0.0136_bc_mutMap_n40_i20,
                    s0.0136_bc_mutMap_n20_i5,
                    s0.0136_bc_mutMap_n20_i10,
                    s0.0136_bc_mutMap_n20_i20,
                    s0.0136_bc_mutMap_n10_i5,
                    s0.0136_bc_mutMap_n10_i10,
                    s0.0136_bc_mutMap_n10_i20)



s0.003_bc_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i5")
s0.003_bc_mutMap_n40_i5$demes <- "40 Demes"
s0.003_bc_mutMap_n40_i5$inds <- "5 Individuals"

s0.003_bc_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i10")
s0.003_bc_mutMap_n40_i10$demes <- "40 Demes"
s0.003_bc_mutMap_n40_i10$inds <- "10 Individuals"

s0.003_bc_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n40_i20")
s0.003_bc_mutMap_n40_i20$demes <- "40 Demes"
s0.003_bc_mutMap_n40_i20$inds <- "20 Individuals"

s0.003_bc_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i5")
s0.003_bc_mutMap_n20_i5$demes <- "20 Demes"
s0.003_bc_mutMap_n20_i5$inds <- "5 Individuals"

s0.003_bc_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i10")
s0.003_bc_mutMap_n20_i10$demes <- "20 Demes"
s0.003_bc_mutMap_n20_i10$inds <- "10 Individuals"

s0.003_bc_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n20_i20")
s0.003_bc_mutMap_n20_i20$demes <- "20 Demes"
s0.003_bc_mutMap_n20_i20$inds <- "20 Individuals"

s0.003_bc_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i5")
s0.003_bc_mutMap_n10_i5$demes <- "10 Demes"
s0.003_bc_mutMap_n10_i5$inds <- "5 Individuals"

s0.003_bc_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i10")
s0.003_bc_mutMap_n10_i10$demes <- "10 Demes"
s0.003_bc_mutMap_n10_i10$inds <- "10 Individuals"

s0.003_bc_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_BC_Map_output_mutMap_n10_i20")
s0.003_bc_mutMap_n10_i20$demes <- "10 Demes"
s0.003_bc_mutMap_n10_i20$inds <- "20 Individuals"

s0.003_bc_mutMap <-rbind( s0.003_bc_mutMap_n40_i5,
                   s0.003_bc_mutMap_n40_i10,
                   s0.003_bc_mutMap_n40_i20,
                   s0.003_bc_mutMap_n20_i5,
                   s0.003_bc_mutMap_n20_i10,
                   s0.003_bc_mutMap_n20_i20,
                   s0.003_bc_mutMap_n10_i5,
                   s0.003_bc_mutMap_n10_i10,
                   s0.003_bc_mutMap_n10_i20)



##################### cline - mutMap

s0.0136_cline_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i5")
s0.0136_cline_mutMap_n40_i5$demes <- "40 Demes"
s0.0136_cline_mutMap_n40_i5$inds <- "5 Individuals"

s0.0136_cline_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i10")
s0.0136_cline_mutMap_n40_i10$demes <- "40 Demes"
s0.0136_cline_mutMap_n40_i10$inds <- "10 Individuals"

s0.0136_cline_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i20")
s0.0136_cline_mutMap_n40_i20$demes <- "40 Demes"
s0.0136_cline_mutMap_n40_i20$inds <- "20 Individuals"

s0.0136_cline_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i5")
s0.0136_cline_mutMap_n20_i5$demes <- "20 Demes"
s0.0136_cline_mutMap_n20_i5$inds <- "5 Individuals"

s0.0136_cline_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i10")
s0.0136_cline_mutMap_n20_i10$demes <- "20 Demes"
s0.0136_cline_mutMap_n20_i10$inds <- "10 Individuals"

s0.0136_cline_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i20")
s0.0136_cline_mutMap_n20_i20$demes <- "20 Demes"
s0.0136_cline_mutMap_n20_i20$inds <- "20 Individuals"

s0.0136_cline_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i5")
s0.0136_cline_mutMap_n10_i5$demes <- "10 Demes"
s0.0136_cline_mutMap_n10_i5$inds <- "5 Individuals"

s0.0136_cline_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i10")
s0.0136_cline_mutMap_n10_i10$demes <- "10 Demes"
s0.0136_cline_mutMap_n10_i10$inds <- "10 Individuals"

s0.0136_cline_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i20")
s0.0136_cline_mutMap_n10_i20$demes <- "10 Demes"
s0.0136_cline_mutMap_n10_i20$inds <- "20 Individuals"

s0.0136_cline_mutMap <-rbind( s0.0136_cline_mutMap_n40_i5,
                        s0.0136_cline_mutMap_n40_i10,
                        s0.0136_cline_mutMap_n40_i20,
                        s0.0136_cline_mutMap_n20_i5,
                        s0.0136_cline_mutMap_n20_i10,
                        s0.0136_cline_mutMap_n20_i20,
                        s0.0136_cline_mutMap_n10_i5,
                        s0.0136_cline_mutMap_n10_i10,
                        s0.0136_cline_mutMap_n10_i20)



s0.003_cline_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i5")
s0.003_cline_mutMap_n40_i5$demes <- "40 Demes"
s0.003_cline_mutMap_n40_i5$inds <- "5 Individuals"

s0.003_cline_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i10")
s0.003_cline_mutMap_n40_i10$demes <- "40 Demes"
s0.003_cline_mutMap_n40_i10$inds <- "10 Individuals"

s0.003_cline_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n40_i20")
s0.003_cline_mutMap_n40_i20$demes <- "40 Demes"
s0.003_cline_mutMap_n40_i20$inds <- "20 Individuals"

s0.003_cline_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i5")
s0.003_cline_mutMap_n20_i5$demes <- "20 Demes"
s0.003_cline_mutMap_n20_i5$inds <- "5 Individuals"

s0.003_cline_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i10")
s0.003_cline_mutMap_n20_i10$demes <- "20 Demes"
s0.003_cline_mutMap_n20_i10$inds <- "10 Individuals"

s0.003_cline_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n20_i20")
s0.003_cline_mutMap_n20_i20$demes <- "20 Demes"
s0.003_cline_mutMap_n20_i20$inds <- "20 Individuals"

s0.003_cline_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i5")
s0.003_cline_mutMap_n10_i5$demes <- "10 Demes"
s0.003_cline_mutMap_n10_i5$inds <- "5 Individuals"

s0.003_cline_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i10")
s0.003_cline_mutMap_n10_i10$demes <- "10 Demes"
s0.003_cline_mutMap_n10_i10$inds <- "10 Individuals"

s0.003_cline_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_cline_output_mutMap_n10_i20")
s0.003_cline_mutMap_n10_i20$demes <- "10 Demes"
s0.003_cline_mutMap_n10_i20$inds <- "20 Individuals"

s0.003_cline_mutMap <-rbind( s0.003_cline_mutMap_n40_i5,
                       s0.003_cline_mutMap_n40_i10,
                       s0.003_cline_mutMap_n40_i20,
                       s0.003_cline_mutMap_n20_i5,
                       s0.003_cline_mutMap_n20_i10,
                       s0.003_cline_mutMap_n20_i20,
                       s0.003_cline_mutMap_n10_i5,
                       s0.003_cline_mutMap_n10_i10,
                       s0.003_cline_mutMap_n10_i20)




##################### trunc - mutMap

s0.0136_trunc_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i5")
s0.0136_trunc_mutMap_n40_i5$demes <- "40 Demes"
s0.0136_trunc_mutMap_n40_i5$inds <- "5 Individuals"

s0.0136_trunc_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i10")
s0.0136_trunc_mutMap_n40_i10$demes <- "40 Demes"
s0.0136_trunc_mutMap_n40_i10$inds <- "10 Individuals"

s0.0136_trunc_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i20")
s0.0136_trunc_mutMap_n40_i20$demes <- "40 Demes"
s0.0136_trunc_mutMap_n40_i20$inds <- "20 Individuals"

s0.0136_trunc_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i5")
s0.0136_trunc_mutMap_n20_i5$demes <- "20 Demes"
s0.0136_trunc_mutMap_n20_i5$inds <- "5 Individuals"

s0.0136_trunc_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i10")
s0.0136_trunc_mutMap_n20_i10$demes <- "20 Demes"
s0.0136_trunc_mutMap_n20_i10$inds <- "10 Individuals"

s0.0136_trunc_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i20")
s0.0136_trunc_mutMap_n20_i20$demes <- "20 Demes"
s0.0136_trunc_mutMap_n20_i20$inds <- "20 Individuals"

s0.0136_trunc_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i5")
s0.0136_trunc_mutMap_n10_i5$demes <- "10 Demes"
s0.0136_trunc_mutMap_n10_i5$inds <- "5 Individuals"

s0.0136_trunc_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i10")
s0.0136_trunc_mutMap_n10_i10$demes <- "10 Demes"
s0.0136_trunc_mutMap_n10_i10$inds <- "10 Individuals"

s0.0136_trunc_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i20")
s0.0136_trunc_mutMap_n10_i20$demes <- "10 Demes"
s0.0136_trunc_mutMap_n10_i20$inds <- "20 Individuals"

s0.0136_trunc_mutMap <-rbind( s0.0136_trunc_mutMap_n40_i5,
                        s0.0136_trunc_mutMap_n40_i10,
                        s0.0136_trunc_mutMap_n40_i20,
                        s0.0136_trunc_mutMap_n20_i5,
                        s0.0136_trunc_mutMap_n20_i10,
                        s0.0136_trunc_mutMap_n20_i20,
                        s0.0136_trunc_mutMap_n10_i5,
                        s0.0136_trunc_mutMap_n10_i10,
                        s0.0136_trunc_mutMap_n10_i20)


s0.003_trunc_mutMap_n40_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i5")
s0.003_trunc_mutMap_n40_i5$demes <- "40 Demes"
s0.003_trunc_mutMap_n40_i5$inds <- "5 Individuals"

s0.003_trunc_mutMap_n40_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i10")
s0.003_trunc_mutMap_n40_i10$demes <- "40 Demes"
s0.003_trunc_mutMap_n40_i10$inds <- "10 Individuals"

s0.003_trunc_mutMap_n40_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n40_i20")
s0.003_trunc_mutMap_n40_i20$demes <- "40 Demes"
s0.003_trunc_mutMap_n40_i20$inds <- "20 Individuals"

s0.003_trunc_mutMap_n20_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i5")
s0.003_trunc_mutMap_n20_i5$demes <- "20 Demes"
s0.003_trunc_mutMap_n20_i5$inds <- "5 Individuals"

s0.003_trunc_mutMap_n20_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i10")
s0.003_trunc_mutMap_n20_i10$demes <- "20 Demes"
s0.003_trunc_mutMap_n20_i10$inds <- "10 Individuals"

s0.003_trunc_mutMap_n20_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n20_i20")
s0.003_trunc_mutMap_n20_i20$demes <- "20 Demes"
s0.003_trunc_mutMap_n20_i20$inds <- "20 Individuals"

s0.003_trunc_mutMap_n10_i5 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i5")
s0.003_trunc_mutMap_n10_i5$demes <- "10 Demes"
s0.003_trunc_mutMap_n10_i5$inds <- "5 Individuals"

s0.003_trunc_mutMap_n10_i10 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i10")
s0.003_trunc_mutMap_n10_i10$demes <- "10 Demes"
s0.003_trunc_mutMap_n10_i10$inds <- "10 Individuals"

s0.003_trunc_mutMap_n10_i20 <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/analysisFiles/AUC_PrecisionRecall_trunc_output_mutMap_n10_i20")
s0.003_trunc_mutMap_n10_i20$demes <- "10 Demes"
s0.003_trunc_mutMap_n10_i20$inds <- "20 Individuals"

s0.003_trunc_mutMap <-rbind( s0.003_trunc_mutMap_n40_i5,
                       s0.003_trunc_mutMap_n40_i10,
                       s0.003_trunc_mutMap_n40_i20,
                       s0.003_trunc_mutMap_n20_i5,
                       s0.003_trunc_mutMap_n20_i10,
                       s0.003_trunc_mutMap_n20_i20,
                       s0.003_trunc_mutMap_n10_i5,
                       s0.003_trunc_mutMap_n10_i10,
                       s0.003_trunc_mutMap_n10_i20)






s0.0136_bc$map <- "BC Map"
s0.0136_bc_mutMap$map <- "BC Map"

s0.003_bc$map <- "BC Map"
s0.003_bc_mutMap$map <- "BC Map"

s0.0136_cline$map <- "Gradient Map"
s0.0136_cline_mutMap$map <- "Gradient Map"

s0.003_cline$map <- "Gradient Map"
s0.003_cline_mutMap$map <- "Gradient Map"

s0.0136_trunc$map <- "Truncated Map"
s0.0136_trunc_mutMap$map <- "Truncated Map"

s0.003_trunc$map <- "Truncated Map"
s0.003_trunc_mutMap$map <- "Truncated Map"

s0.0136 <- rbind( s0.0136_bc,
                  s0.0136_cline,
                  s0.0136_trunc)
s0.0136$selection <- "s = 0.0136"
s0.0136$mutation <- "Constant Mutation Rate"

s0.003 <- rbind( s0.003_bc,
                 s0.003_cline,
                 s0.003_trunc)
s0.003$selection <- "s = 0.003"
s0.003$mutation <- "Constant Mutation Rate"


s0.0136_mutMap <- rbind( s0.0136_bc_mutMap,
                  s0.0136_cline_mutMap,
                  s0.0136_trunc_mutMap)
s0.0136_mutMap$selection <- "s = 0.0136"
s0.0136_mutMap$mutation <- "Variable Mutation Rate"

s0.003_mutMap <- rbind( s0.003_bc_mutMap,
                 s0.003_cline_mutMap,
                 s0.003_trunc_mutMap)
s0.003_mutMap$selection <- "s = 0.003"
s0.003_mutMap$mutation <- "Variable Mutation Rate"


all_mutMap <- rbind(s0.003_mutMap, 
             s0.0136_mutMap)

all_uniform <- rbind(s0.003, 
                    s0.0136)

all <- rbind(all_uniform, all_mutMap)

stats <- c("geno_k_tau_p_value",
           "LEA",
           "RDA",
           "BF(dB)",
           "top_candidate_test" ,
           "Z_geno_k_tau_p_value", 
           "Z_LEA",
           "Z_RDA",
           "Z_BF(dB)")


WZA_stats <- c("Z_geno_k_tau_p_value", 
           "Z_LEA",
           "Z_RDA",
           "Z_BF(dB)")

labels_main_figs <- c(expression("Kendall's "*tau),
                      expression("LFMM - LEA"),
                      expression("RDA"),
                      expression("BayPass"),
                      expression("Top-Candidate"),
                      expression(bold("WZA"[tau])),
                      expression("WZA"["LFMM"]),
                      expression("WZA"["RDA"]),
                      expression("WZA"["BayPass"]))

main_stats <- c("geno_k_tau_p_value",
                  "LEA",
                  "RDA",
                  "BF(dB)",
                  "top_candidate_test" ,
                  "Z_geno_k_tau_p_value")

all$statistic  <- factor(all$stat,
                         levels = stats,
                         labels = labels_main_figs)

all$individuals  <- factor(all$inds,
                         levels = c("5 Individuals",
                                    "10 Individuals",
                                    "20 Individuals"))

parse.labels <- function(x) parse(text = x)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "white")

str(all)

figure_3_data_constant_weak <- all[(all$inds == "20 Individuals")&
                                     (all$mutation == "Constant Mutation Rate")&
                                     (all$stat%in%main_stats)&
                                   (all$selection=="s = 0.003"),]
figure_3_data_constant_strong <- all[(all$inds == "20 Individuals")&
                                     (all$mutation == "Constant Mutation Rate")&
                                     (all$stat%in%main_stats)&
                                   (all$selection=="s = 0.0136"),]
figure_3_data_variable_weak <- all[(all$inds == "20 Individuals")&
                                       (all$mutation == "Variable Mutation Rate")&
                                       (all$stat%in%main_stats)&
                                       (all$selection=="s = 0.003"),]
figure_3_data_variable_strong <- all[(all$inds == "20 Individuals")&
                                (all$mutation == "Variable Mutation Rate")&
                                (all$stat%in%main_stats)&
                                (all$selection=="s = 0.0136"),]


figure_3B_variable_strong <- ggplot(data = figure_3_data_variable_strong, 
       aes(x = statistic,
           y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

figure_3B_variable_weak <- ggplot(data = figure_3_data_variable_weak, 
                    aes(x = statistic,
                        y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

figure_3B_constant_strong <- ggplot(data = figure_3_data_constant_strong, 
                                  aes(x = statistic,
                                      y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  scale_color_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))


figure_3B_constant_weak <- ggplot(data = figure_3_data_constant_weak, 
                                  aes(x = statistic,
                                      y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
    stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))


ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Figure_3_variable_weak.pdf",
       plot = figure_3B_variable_weak,
       width = 9.08,
       height = 5.71)
ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Figure_3_variable_strong.pdf",
       plot = figure_3B_variable_strong,
       width = 9.08,
       height = 5.71)
ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Figure_3_constant_weak.pdf",
       plot = figure_3B_constant_weak,
       width = 9.08,
       height = 5.71)
ggsave("~/work/GEA/simulations/SimulationStudy/Plots/Figure_3_constant_strong.pdf",
       plot = figure_3B_constant_strong,
       width = 9.08,
       height = 5.71)



suppfigure_sampleSize_data_constant_weak <- all[(all$mutation == "Constant Mutation Rate")&
                                     (all$stat%in%main_stats)&
                                     (all$selection=="s = 0.003"),]
suppfigure_sampleSize_data_constant_strong <- all[(all$mutation == "Constant Mutation Rate")&
                                       (all$stat%in%main_stats)&
                                       (all$selection=="s = 0.0136"),]


suppfigure_sampleSize_constant_weak <- ggplot(data = suppfigure_sampleSize_data_constant_weak, 
                                  aes(x = statistic,
                                      y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map*individuals~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =10),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

suppfigure_sampleSize_constant_weak <- ggplot(data = suppfigure_sampleSize_data_constant_strong, 
                                              aes(x = statistic,
                                                  y = auc_pr))+
  annotate(geom = "rect", xmin = 5.5 , xmax = 6.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map*individuals~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =10),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("~/work/GEA/simulations/SimulationStudy/Plots/suppfigure_sampleSize_constant_strong.pdf",
       plot = suppfigure_sampleSize_constant_weak,
       width = 9.08,
       height = 5.71*2)


ggsave("~/work/GEA/simulations/SimulationStudy/Plots/suppfigure_sampleSize_constant_weak.pdf",
       plot = suppfigure_sampleSize_constant_weak,
       width = 9.08,
       height = 5.71*2)




#### Neutral histograms

wza_bc_neutral_data_var <- list()
wza_cline_neutral_data_var <- list()
wza_trunc_neutral_data_var <- list()

wza_bc_neutral_data_con <- list()
wza_cline_neutral_data_con <- list()
wza_trunc_neutral_data_con <- list()

for (i in 1:30){
  wza_tmp_bc_var <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.003/BC_Map_output_mutMap_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))
  wza_tmp_cline_var <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.003/cline_output_mutMap_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))
  wza_tmp_trunc_var <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.003/trunc_output_mutMap_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))

  wza_bc_neutral_data_var[[i]] <- wza_tmp_bc_var[wza_tmp_bc_var$gene%in%paste0("gene",800:1000),]
  wza_cline_neutral_data_var[[i]] <- wza_tmp_cline_var[wza_tmp_cline_var$gene%in%paste0("gene",800:1000),]
  wza_trunc_neutral_data_var[[i]] <- wza_tmp_trunc_var[wza_tmp_trunc_var$gene%in%paste0("gene",800:1000),]

  wza_tmp_bc_con <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.0136/BC_Map_output_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))
  wza_tmp_cline_con <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.0136/cline_output_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))
  wza_tmp_trunc_con <- read.csv(paste0("~/work/GEA/simulations/SimulationStudy/s0.0136/trunc_output_n40_i20/",i,"_WZA_geno_k_tau_p_value.csv"))
  
  wza_bc_neutral_data_con[[i]] <- wza_tmp_bc_con[wza_tmp_bc_con$gene%in%paste0("gene",800:1000),]
  wza_cline_neutral_data_con[[i]] <- wza_tmp_cline_con[wza_tmp_cline_con$gene%in%paste0("gene",800:1000),]
  wza_trunc_neutral_data_con[[i]] <- wza_tmp_trunc_con[wza_tmp_trunc_con$gene%in%paste0("gene",800:1000),]
  
    #  snp_neutral_data[[i]] <- snp_tmp[snp_tmp$gene%in%paste0("gene",800:1000),]
}


wza_bc_neutral_var <- do.call(rbind, wza_bc_neutral_data_var)
wza_bc_neutral_var$map <- "BC Map"
wza_cline_neutral_var <- do.call(rbind, wza_cline_neutral_data_var)
wza_cline_neutral_var$map <- "Gradient Map"
wza_trunc_neutral_var <- do.call(rbind, wza_trunc_neutral_data_var)
wza_trunc_neutral_var$map <- "Truncated Map"

wza_neutral_var <- rbind(wza_bc_neutral_var,
                         wza_cline_neutral_var,
                         wza_trunc_neutral_var)
wza_neutral_var$mutRate <- "Variable Mutation Rate"

wza_bc_neutral_con <- do.call(rbind, wza_bc_neutral_data_con)
wza_bc_neutral_con$map <- "BC Map"
wza_cline_neutral_con <- do.call(rbind, wza_cline_neutral_data_con)
wza_cline_neutral_con$map <- "Gradient Map"
wza_trunc_neutral_con <- do.call(rbind, wza_trunc_neutral_data_con)
wza_trunc_neutral_con$map <- "Truncated Map"

wza_neutral_con <- rbind(wza_bc_neutral_con,
                         wza_cline_neutral_con,
                         wza_trunc_neutral_con)
wza_neutral_con$mutRate <- "Constant Mutation Rate"

wza_neutral <- rbind(wza_neutral_con, wza_neutral_var)

neutral_histograms <- ggplot(data = wza_neutral,
       aes(x = Z_pVal))+
  geom_histogram(fill = "white", 
                 col = "black",
                 boundary = 1)+
  scale_x_continuous(limits = c(0,1))+
  facet_grid(map ~ mutRate)+
  guides(fill = "none",
         color = "none")+
  ylab("Count")+
  xlab(expression(italic(p)*"-value"))+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8)+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8)+
  theme_bw()+
  
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title = element_text(size =13))


ggsave("~/work/GEA/simulations/SimulationStudy/Plots/neutralHistograms.pdf",
       plot = neutral_histograms,
       width = 9.08,
       height = 5.71)


## Different flavours of the WZA


suppfigure_WZAstat_constant_weak <- all[(all$inds == "20 Individuals")&
                                          (all$mutation == "Constant Mutation Rate")&
                                          (all$stat%in%WZA_stats)&
                                          (all$selection=="s = 0.003"),]
suppfigure_WZAstat_constant_strong <- all[(all$inds == "20 Individuals")&
                                            (all$mutation == "Constant Mutation Rate")&
                                            (all$stat%in%WZA_stats)&
                                            (all$selection=="s = 0.0136"),]



suppfigure_WZAstat_constant_strong_plot  <- ggplot(data = suppfigure_WZAstat_constant_strong, 
                                                   aes(x = statistic,
                                                       y = auc_pr))+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8[6:9])+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8[6:9])+
  theme_bw()+
  
  theme(axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))



suppfigure_WZAstat_constant_weak_plot  <- ggplot(data = suppfigure_WZAstat_constant_weak, 
                                                   aes(x = statistic,
                                                       y = auc_pr))+
  geom_violin(aes(fill = statistic), alpha = 0.5)+
  geom_jitter(aes(fill = statistic), pch =21,col = "black", height = 0, width = 0.1)+
  stat_summary(fun="mean", geom="segment", 
               mapping=aes(xend=..x.. + 0.5, yend=..y.., col = statistic), lwd = 1)+
  facet_grid(map~demes)+
  coord_cartesian(ylim = c(0,1))+
  scale_x_discrete(labels = parse.labels)+
  guides(fill = "none",
         color = "none")+
  xlab(NULL)+
  ylab("Area Under Precision Recall Curve")+
  scale_color_manual(labels = parse.labels,
                     values =colorBlindGrey8[6:9])+
  scale_fill_manual(labels = parse.labels,
                    values =colorBlindGrey8[6:9])+
  theme_bw()+
  
  theme(axis.text.x = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =13),
        axis.title.y = element_text(face = "bold", 
                                    size =13))


ggsave("~/work/GEA/simulations/SimulationStudy/Plots/suppfigure_WZAstat_constant_strong.pdf",
       plot = suppfigure_WZAstat_constant_strong_plot,
       width = 9.08,
       height = 5.71)

ggsave("~/work/GEA/simulations/SimulationStudy/Plots/suppfigure_WZAstat_constant_weak.pdf",
       plot = suppfigure_WZAstat_constant_weak_plot,
       width = 9.08,
       height = 5.71)


# Power plots

s0.0136_power <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.0136/powerAnalysis.csv")
s0.0136_power[s0.0136_power$method=="RDA",]$fdr = 0
s0.0136_power <- s0.0136_power[s0.0136_power$mutation=="variation",]
s0.0136_power_melted <- melt(s0.0136_power, id = c("method","map","mutation","rep","demes","inds"))
s0.0136_power_melted$deme_id <- paste0(s0.0136_power_melted$demes, " Demes")
s0.0136_power_melted$map <-factor(s0.0136_power_melted$map,
                                 levels = c("BC","cline","trunc"),
                                 labels = c("BC Map",
                                            "Gradient Map",
                                            "Truncated Map"))



main_stats <- c("LEA","RDA","BF","WZA")
labels_for_figs <- c(expression("LFMM - LEA"),
                     expression("RDA"),
                     expression("BayPass"),
                     expression(bold("WZA"[tau])))

s0.0136_power_melted$statistic  <- factor(s0.0136_power_melted$method,
                                        levels = main_stats,
                                        labels = labels_for_figs)
s0.0136_power_melted$variable <- factor(s0.0136_power_melted$variable,
                                       levels = c("fdr","power"),
                                       labels = c("FDR","Power"))


parse.labels <- function(x) parse(text = x)

colorBlindGrey8   <- c("#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7", "white")


s0.0136_powerPlot <- ggplot(data = s0.0136_power_melted[s0.0136_power_melted$inds==20,], aes(x = statistic, shape = variable, y = value, fill = statistic))+
  annotate(geom = "rect", xmin = 3.5 , xmax = 4.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  facet_grid(map~deme_id)+
  geom_point(position = position_dodge(width = 1), size = 1.4)+
  stat_summary(aes(shape = variable),
               position = position_dodge(width = 1),
               size = 1)+
  annotate("text", x = 2-.25, y = .075, label = "*", size = 10)+
  scale_fill_manual(values =colorBlindGrey8)+
  scale_shape_manual("",values =c(21,24))+
  coord_cartesian(ylim = c(0,1))+
  guides(fill = "none")+
  scale_x_discrete(labels = parse.labels)+
  theme_bw()+
  xlab(NULL)+
  ylab(NULL)+
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =10),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("/home/booker/work/GEA/simulations/SimulationStudy/Plots/suppfigure_power_plot_strong.pdf",
       s0.0136_powerPlot,
       width = 9.08,
       height = 5.71 )

s0.003_power <- read.csv("~/work/GEA/simulations/SimulationStudy/s0.003/powerAnalysis.csv")
s0.003_power[s0.003_power$method=="RDA",]$fdr = 0
s0.003_power <- s0.003_power[s0.003_power$mutation=="variation",]
s0.003_power_melted <- melt(s0.003_power, id = c("method","map","mutation","rep","demes","inds"))
s0.003_power_melted$deme_id <- paste0(s0.003_power_melted$demes, " Demes")
s0.003_power_melted$map <-factor(s0.003_power_melted$map,
                                 levels = c("BC","cline","trunc"),
                                 labels = c("BC Map",
                                            "Gradient Map",
                                            "Truncated Map"))



main_stats <- c("LEA","RDA","BF","WZA")
labels_for_figs <- c(expression("LFMM - LEA"),
                     expression("RDA"),
                     expression("BayPass"),
                     expression(bold("WZA"[tau])))


s0.003_power_melted$statistic  <- factor(s0.003_power_melted$method,
                                         levels = main_stats,
                                         labels = labels_for_figs)
s0.003_power_melted$variable <- factor(s0.003_power_melted$variable,
                                       levels = c("fdr","power"),
                                       labels = c("FDR","Power"))

s0.003_powerPlot <- ggplot(data = s0.003_power_melted[s0.003_power_melted$inds==20,], aes(x = statistic, shape = variable, y = value, fill = statistic))+
  annotate(geom = "rect", xmin = 3.5 , xmax = 4.5 , ymin = -1.1 , ymax = 40, fill = "grey", alpha = 0.5)+
  facet_grid(map~deme_id)+
  geom_point(position = position_dodge(width = 1), size = 1.4)+
  stat_summary(aes(shape = variable),
               position = position_dodge(width = 1),
               size = 1)+
  annotate("text", x = 2-.25, y = .075, label = "*", size = 10)+
  scale_fill_manual(values =colorBlindGrey8)+
  scale_shape_manual("",values =c(21,24))+
  coord_cartesian(ylim = c(0,1))+
  guides(fill = "none")+
  scale_x_discrete(labels = parse.labels)+
  theme_bw()+
  xlab(NULL)+
  ylab(NULL)+
  theme(axis.text.x = element_text(angle = 80, 
                                   hjust = 1, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size =10),
        axis.title.y = element_text(face = "bold", 
                                    size =13))

ggsave("/home/booker/work/GEA/simulations/SimulationStudy/Plots/suppfigure_power_plot_weak.pdf",
       s0.003_powerPlot,
       width = 9.08,
       height = 5.71 )



