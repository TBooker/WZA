import numpy as np
import pandas as pd
import gzip
import argparse
import scipy.stats

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="A script to compute GEA using Kendall's tau from BayPass input files")
	parser.add_argument("--geno",
			required = True,
			dest = "geno",
			type = str,
			help = "The genotype data in BayPass format - see the manual of that software for info")
	parser.add_argument("--env",
			required = True,
			dest = "env",
			type = str,
			help = "The file of environmental variables (assumed to be in the same order as the genotype file)")
	parser.add_argument("--meta",
			required = True,
			dest = "meta",
			type = str,
			help = "The metadata of the SNPs that you are analysing (see vignette for details)")

	parser.add_argument("--output",
			required = True,
			dest = "output",
			type = str,
			help = "What name do you want to give the output file?")

	args = parser.parse_args()

# First thing we'll do is to read in the environmental data
# this next line is a bit ugly, but it gets the job done.
# It reads in the environmental data, which is stored on a single line
# Then does some simple processing to convert it from text to a numerical vector

	env_vector =  np.array( list( map(float, open(args.env).readlines()[0].strip().split(" ") ) ))


	metadata = pd.read_csv(args.meta)



	count = 0

	results = []

## Open up the genotype file and iterate over all SNPs in the dataset
	with gzip.open(args.geno) as file:
		for line in file:
			count +=1
## The following line is quite similar to the one above
# It includes ".decode("utf-8")" which converts the text in byte format
# into a string characters
			raw_snp_freqs = np.array([float(i) for i in line.decode("utf-8").strip("\n").split(" ") if i!=""])

## Now calculate the allele frequnecy for each population
			snp_freqs = np.array([raw_snp_freqs[k]/(raw_snp_freqs[k]+raw_snp_freqs[k+1]) for k in range(0,raw_snp_freqs.shape[0],2)])

## Now calculate the mean allele frequncy across populations
			p_bar = snp_freqs.mean()
			if p_bar > 0.5:
				maf = 1 - p_bar
			else:
				maf = p_bar

			geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(env_vector, snp_freqs)

			outline = [geno_k_tau, geno_k_tau_p_value, maf]
			results.append(outline)


	if len(results) != metadata.shape[0]:
		print("metadata and the genotype data are different shapes")
		return
	geno_output = pd.DataFrame(results, columns = ["K_tau","K_tau_p", "MAF"])

	all_output = pd.concat([geno_output,metadata], axis = 1)

	all_output.to_csv(args.output,
						index = False,
						sep = "\t")

if __name__ == "__main__":
	main()
