import pandas as pd
import scipy.stats
import numpy as np
import argparse, tskit
import matplotlib.pyplot as pyplot

## This script will take the output from many GEA simulation runs
## combine them together and perform multiple tests

def main():
	
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--GEA_dir","-g", 
			required = True,
			dest = "GEA_dir",
			type = str, 
			help = "The directory with the GEA results (saved as CSV files)")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "What name do you want to give to the output file?")
	args = parser.parse_args()
