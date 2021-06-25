import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse

## Script to convert the Pine Bayes Factor results into empirical p-values
def main():

## Define command line args

	parser = argparse.ArgumentParser(description="")

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

			help = "If you want to analyse just a single environment, give it here [DD_0]",
			
			default = "DD_0")
			
	parser.add_argument("--bay", 

			required = False,

			dest = "bay",

			action = "store_true", 

			help = "[OPTIONAL] Give this flag if the analysis files are BayesFactors from bayEnv.")
			
			
	
main()


