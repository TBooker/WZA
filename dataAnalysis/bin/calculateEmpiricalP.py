import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse

def main():

## Define command line args

	parser = argparse.ArgumentParser(description="A quick script to calculate empirical p-values")

	parser.add_argument("--correlations", "-c", 

			required = True,

			dest = "correlations",

			type = str, 

			help = "The file containing the correlations")

	parser.add_argument("--output", 

			required = True,

			dest = "output",

			type = str, 

			help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

	parser.add_argument("--env", 

			required = False,

			dest = "env",

			type = str, 

			help = "[OPTIONAL] If you want to analyse just a single environment, give it here. The default is to run the analysis on all CoAdapTree envs")

	parser.add_argument("--maf", 

			required = False,

			dest = "maf",

			type = float, 

			help = "[OPTIONAL] If you want to apply a MAF filter specify that here",
			default = 0)
			
	args = parser.parse_args()
	

## For all environmental variables:
## MAT MWMT MCMT TD MAP MSP AHM SHM DD_0 DDS NFFD bFFP FFP PAS EMT EXT Eref CMD
	if args.env:	
		envs = [args.env]
	else:
		envs =["MAT","MWMT","MCMT","TD",
			"MAP","MSP","AHM","SHM",
			"DD_0","DD5","NFFD","bFFP",
			"FFP","PAS","EMT","EXT","Eref","CMD"]
## Make sure the env you gave is actually a in the dataset
	if args.env and args.env not in envs:
		print(args.env + " is not a valid environment")
		return
	
### This is going to use a LOT of memory and might not be feasible in some cases
	big_CSV = pd.read_csv( args.correlations , sep  = "\t")

	for environment in envs:
#	for environment in ["AHM"]:
	## Extract the SNPs corresponding to the environment in question
		temp = big_CSV[big_CSV["env"] == environment].copy()
	## Apply a MAF 
		if args.maf >0:
			temp = temp[temp.MAF > args.maf]
	
	## Calculate the empirical p-values
		temp["empirical_pvalue"] = temp["pvalue"].rank()/temp.shape[0]
	## If there were no SNPs corresponding to the current env then don't write a file
		if temp.shape[0] == 0:
			continue
		print(temp.shape)
		temp.to_csv( environment + "_empiricalP_" + args.correlations, index = False, sep = '\t')

main()
