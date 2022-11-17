## A script to subsample SNPs from a window and to use those as input to the WZA...

import pandas as pd
import argparse



def main():

## Define command line args
	parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

	parser.add_argument("--correlations", "-c",
			required = True,
			dest = "correlations",
			type = str,
			help = "The file containing the correlations")

	parser.add_argument("--summary_stat", "-s",
			required = True,
			dest = "summary_stat",
			type = str,
			help = "The name of the column you are analysing")

	parser.add_argument("--window", "-w",
			required = True,
			dest = "window",
			type = str,
			help = "The name of column containing the windows you want to analyse")

	parser.add_argument("--output",
			required = True,
			dest = "output",
			type = str,
			help = "The name of the output file")

	parser.add_argument("--sample_snps",
			required = False,
			dest = "sample_snps",
			type = int,
			help = "[OPTIONAL] Give the number of SNPs you want to downsample to. Give -1 if you want to use the 75th percentile of the number of SNPs. Note that calculating the median within the script is slow, so you may want to run a dummy analysis, get the median number of SNPs then use that explicitly.",
			default = 0)

	parser.add_argument("--resamples",
			required = False,
			dest = "resamples",
			type = int,
			help = "[OPTIONAL] Give the number of times you want to resample WZA scores when the number of SNPs exceeds the sample_snps threshold. [100]",
			default = 100)

	parser.add_argument("--min_snps",
			required = False,
			dest = "min_snps",
			type = int,
			help = "[OPTIONAL] Give the minimum number of SNPs you are willing to analyse per window. Default is 3",
			default = 3)

	parser.add_argument("--large_i_small_p",
			required = False,
			action = "store_true",
			help = "[OPTIONAL] Give this flag if extreme values of the summary stat you're using are large values.")

	parser.add_argument("--top_candidate_threshold",
			required = False,
			dest = "top_candidate_threshold",
			type = float,
			help = "[OPTIONAL] Give the percentile threshold you want to use for the top-candidate test (from Yeaman et al 2016 - Science)",
			default = 99)

	parser.add_argument("--verbose", "-v",
			required = False,
			action = "store_true",
			help = "[OPTIONAL] Give this flag if you want to run the script in verbose mode")

	parser.add_argument("--MAF",
			required = False,
			dest = "MAF",
			type = str,
			help = "[OPTIONAL] Give the name of the MAF colummn - if it's not explicitly given in the input dataframe")

	parser.add_argument("--sep",
			required = False,
			dest = "sep",
			type = str,
			default = "\t",
			help = "What separator do you use in your file?")

	parser.add_argument("--retain",
			required = False,
			dest = "retain",
			nargs="+",
			type = str,
			help = "Give any columns that you want to add to the output file")

	args = parser.parse_args()

	csv = pd.read_csv(args.correlations, sep = args.sep)

	print(csv)



main()
