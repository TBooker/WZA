import msprime, pyslim, random, allel, glob
import pandas as pd
import scipy.stats
from scipy.spatial import distance_matrix
import numpy as np
import argparse, tskit
import matplotlib.pyplot as plt



def summary_stats(tree_sequence_file):
#	ts = pyslim.load(tree_sequence_file)
	ts = pyslim.load(tree_sequence_file)
	mutated_tree = msprime.mutate(ts, 1e-7)
#	muts = len( [ v for v  in mutated_tree.variants() ] )
# Get the genotype matrix, ready for using sci-kit.allel
	msprime_genotype_matrix = mutated_tree.genotype_matrix()
# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
	haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )

	genotype_array = haplotype_array.to_genotypes(ploidy=2)

## Calculate Diversity
	pi = mutated_tree.diversity(windows =[0,1e6,1e6+3000])
#	print(pi, genotype_array.shape)
	subpopulations = [ [y for y in range(x, x+100)] for x in range(0,genotype_array.shape[1],100)]
	
#	print(len(individuals), genotype_array.shape)
	subpops = np.array(subpopulations)[np.random.choice(len(subpopulations),10, replace = False)]
	mean_fst = allel.average_weir_cockerham_fst(genotype_array, blen = 10000, subpops=subpops)
#	print(mean_fst)
	return(pi[0], mean_fst[0])

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
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
	parser.add_argument("--noTime", 
			dest = "noTime",
			action = "store_true", 
			help = "Use this flag if want to ignore the time stuff")
	args = parser.parse_args()

	output = open(args.output, "w")
	output.write(','.join([ "rep", "gen", "pi", "Fst", "correlation"]) + "\n")
	
	count = 0
	for i in glob.glob(args.trees + "/*trees"):
		count +=1
		print(i)
		name = i.split("/")[-1].split(".")
		direc = "/".join(i.split("/")[0:-1])
		rep = name[0]


		if args.noTime:
			gen = "XX"
			phen = pd.read_csv(i.split("trees")[0] + "phen.txt")
		else:
			gen = name[1]
			phen = pd.read_csv(direc + "/" + gen +"." + rep + ".checkMeOutChump.phen.txt")
		if gen == "200102":
			gen = "210102"

		pi, Fst = summary_stats(i)
# I made a dumb error in the SLiM code that caused me to name the phenotype feils for generation 210102 incorrectely.
		
		phen_env_corr = scipy.stats.kendalltau(phen.phen, phen.opt)
		cor = phen_env_corr.correlation
		print( [rep, gen, str(pi), str(Fst), str(cor)] )
		output.write(",".join([rep, gen, str(pi), str(Fst), str(cor)]) + "\n")
#		if count == 5:break
	output.close()

main()
