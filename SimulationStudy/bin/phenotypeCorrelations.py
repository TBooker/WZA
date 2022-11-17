## A short little script to get the phen, env correlation across a bunch of simulations
import glob
import pandas as pd
import scipy.stats
import numpy as np
import argparse

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--phen", 
			required = True,
			dest = "phen",
			type = str, 
			help = "Directory with the phenotype files")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "Give the output files name")
	args = parser.parse_args()
	
	corDat = []
	count = 0

	parameter_combinations = ["0.3_0.001_100", "0.3_0.0001_100", "1_0.001_100","1_0.0001_100", "0.3_0.001_200", "0.3_0.0001_200", "1_0.001_200","1_0.0001_200"]
	results = []
	for i in parameter_combinations:
		print(i)
		sig_a = i.split("_")[0]
		Ua = i.split("_")[1]
		Vs = i.split("_")[2]
		for  j in glob.glob(args.phen + "*"+ i + ".phen.txt"):
			rep = j.split("/")[-1].split("_")[0]
			phenOpt = pd.read_csv(j)
			pop_k_tau, pop_k_tau_p_value = scipy.stats.kendalltau(phenOpt.opt, phenOpt.phen)
			results.append( [ rep, sig_a, Ua, Vs, pop_k_tau, pop_k_tau_p_value, i ] )
	output = pd.DataFrame(results,
				columns = ["rep","sig_a","Ua","Vs","tau","pVal","params"] )
	output.to_csv(args.output, index = False) 

"""
	for (i in part1){
	  count= count+1
	  data_vec = rep(0,100)
	  for (j in 1:98){ 
		cat(paste("/media/booker/HOWDY/GEA/G.0.8_selected2/",j,"_",i,".phen.txt",sep = ""))
		x <- read.csv(paste("/media/booker/HOWDY/GEA/G.0.8_selected2/",j,"_",i,".phen.txt",sep = ""))
		data_vec[j] = cor(x$phen, x$opt)
	  }
	  nameList = strsplit(i,"\\_")[[1]]
	  corDat[[count]]  = data.frame(sig_a=nameList[1],U_a=nameList[2], Vs=nameList[3], rep = 1:100, cor = data_vec)
	}

	correlations <- do.call(rbind, corDat)

	correlations$sig_a <- factor(as.factor(correlations$sig_a),
								 levels = c("0.3",
											"0.5",
											"1"), 
								 labels = c(expression(italic(sigma[a]^"2")*" = 0.3"),
											expression(italic(sigma[a]^"2")*" = 0.5"),
											expression(italic(sigma[a]^"2")*" = 1.0" )
								 ))

	correlations$U_a <- factor(as.factor(correlations$U_a),
							   levels = c("0.001",
										  "0.0001"), 
							   labels = c(expression(italic(U[a])*" = 0.001"),
										  expression(italic(U[a])*" = 0.0001" )
							   ))

	correlations$Vs <- factor(as.factor(correlations$Vs),
							  levels = c("100",
										 "200"), 
							  labels = c(expression(V[s]*" = 100"),
										 expression(V[s]*" = 200" )
							  ))
"""

main()
