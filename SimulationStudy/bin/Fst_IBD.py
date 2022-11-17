import msprime, pyslim, random, allel, glob, sys
import pandas as pd
import scipy.stats
from scipy.spatial import distance_matrix
import numpy as np
import argparse, tskit
import matplotlib.pyplot as plt




def Fst_IBD(trees):
	ts = pyslim.load(trees)
	mutated_tree = msprime.mutate(ts, 1e-8, keep = False)
#	muts = len( [ v for v  in mutated_tree.variants() ] )
# Get the genotype matrix, ready for using sci-kit.allel
	msprime_genotype_matrix = mutated_tree.genotype_matrix()
# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
	haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )

	genotype_array = haplotype_array.to_genotypes(ploidy=2)
	print(genotype_array.shape)
## Calculate Diversity
#	pi = mutated_tree.diversity(windows =[0,1e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,10e6,10e6+1])
## Calculate Tajima's D
#	ac = genotype_array.count_alleles()
#	TD = allel.tajima_d(ac)
#	print(TD)


	row = np.random.choice(13)
	pairs = [[row, row+(14*i)]  for i in range(14)]

	subpopulations = [ [y for y in range(x, x+100)] for x in range(0,genotype_array.shape[1],100)]

	subpops = np.array(subpopulations)[np.random.choice(len(subpopulations),10, replace = False)]

	mean_fst = allel.average_weir_cockerham_fst(genotype_array, blen = 10000, subpops= subpops)

	rep = trees.split("/")[-1].split("_")[0]

	output = []

	output.append( [ str(rep), str(int(-1)), str(mean_fst[0]) ])

	for p in pairs:
		dist = (p[1] - p[0])/14
		if dist == 0: continue
		subpops = np.array(subpopulations)[p]

		mean_fst = allel.average_weir_cockerham_fst(genotype_array, blen = 1000, subpops=subpops)
		output.append( [ str(rep), str(int(dist)), str(mean_fst[0]) ])
		print(p, mean_fst)
#		output.write( ",".join( str(rep), str(int(dist)), str(mean_fst[0]) ) + "\n")
	return(output)

def main():
	parser = argparse.ArgumentParser(description="Let's look at Fst between differing pairs of demes in the 2D stepping-stone")
	parser.add_argument("--trees",
			required = True,
			dest = "trees",
			type = str,
			help = "Directory with the coalescent trees for the simulation")
	parser.add_argument("--output",
			required = True,
			dest = "output",
			type = str,
			help = "What name do you want to give to the output file? [DON'T GIVE FILE EXTENSION, THE SCRIPT MAKES TWO OUTPUT FILES")
	args = parser.parse_args()
	output_data = Fst_IBD(args.trees)
	pd.DataFrame(output_data, columns = ["rep","dist","Fst"]).to_csv(args.output, index = False)

main()
