import msprime, pyslim, random, allel, glob, itertools
import pandas as pd
import scipy.stats
from scipy.spatial import distance_matrix
import numpy as np
import argparse, tskit
import matplotlib.pyplot as plt




# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
# For item i in a range that is a list of length l,
	for i in range(0, len(l), n):
		# Create an index range for l of n items:
		yield l[i:i+n]

def getLD(tree_sequence, output_name):
	
	pops = {} 
	for i in tree_sequence.individuals_alive_at(0):
		if tree_sequence.individual(i).population == 999: continue

		try: 
			pops[tree_sequence.individual(i).population].append( i )
		except KeyError:
			pops[tree_sequence.individual(i).population] = [i]
	
# Let's take a sample of 10 individuals from 1 pops and make a new tree from them
	sampled_indivs_1 = []
	for i in np.random.choice( list( pops.keys() ), 1, replace = False):
		for j in np.random.choice( np.array(pops[i]), 10 , replace = False):
			sampled_indivs_1.extend( tree_sequence.individual(j).nodes )
	
	diploids_per_pop = 2
# Let's take a sample of 10 individuals from 20 pops and make a new tree from them
	sampled_indivs_2 = []
	for i in np.random.choice( list( pops.keys() ), 20, replace = False):
		for j in np.random.choice( np.array(pops[i]), 2 , replace = False):
			sampled_indivs_2.extend( tree_sequence.individual(j).nodes )
		
	
	
	r2_tree_1 = tree_sequence.simplify(sampled_indivs_1)
	mut_tree_1 = msprime.mutate(r2_tree_1, 1e-7)

	r2_tree_2 = tree_sequence.simplify(sampled_indivs_2)
	mut_tree_2 = msprime.mutate(r2_tree_2, 1e-7)


#	print(ts.sample_size)
## Sprinkle mutations onto the coalescent tree
#	print('Sprinkling mutations onto trees')

#	print("sprinkle 1")
#	sprinkled = msprime.mutate(r2_tree2, rate= mut_rate, keep=True)

	LD_intervals = [[9000000,9010000-1], [9200000,9210000-1], [9400000,9410000-1], [9600000,9610000-1], [9800000,9810000-1]]


	LD_intervals = [[7260000,7270000-1]]

	LD_output_within = open("within_"+output_name, "w")
	for interval in LD_intervals:
		positions = [i.position for i in mut_tree_1.variants() if i.position >= interval[0] and i.position <= interval[1]]
		LD_interval = mut_tree_1.keep_intervals(np.array([ interval ])) 
		LD = tskit.LdCalculator(LD_interval)
		LD_mat = LD.r2_matrix()
		freq_dict = {}
		for variant in LD_interval.variants():
## Calculating the genotype freqs this way assumes that the ref and alt are coded as 1s and 0s, repectively. I'll keep that for now, but will change in a bit
			allele_freq = ((variant.genotypes[::2] + variant.genotypes[1::2])/2).mean()
			freq_dict[variant.position] = allele_freq

		print(len(freq_dict.keys()))
		print(len(positions))
		
		for i in range(len(positions)):
			for j in range(len(positions)):
				if j>=i:
					continue
				else:
					if freq_dict[positions[i]] < 0.1 or freq_dict[positions[j]] < 0.1: continue 
					pw_distances_ij = abs( positions[i] - positions[j] )
					if pw_distances_ij > 10000: continue
					LD_output_within.write(",".join([str(int(pw_distances_ij)), str(LD_mat[i,j]) ] ) + "\n")
					
	LD_output_within.close()
	print("Now getting 'LD' between pops")
	LD_output_between = open("between_"+output_name, "w")

	for interval in LD_intervals:
		positions = [i.position for i in mut_tree_2.variants() if i.position >= interval[0] and i.position <= interval[1]]
		LD_interval = mut_tree_2.keep_intervals(np.array([ interval ])) 

		freq_dict = {}
		
		for variant in LD_interval.variants():

## Calculating the genotype freqs this way assumes that the ref and alt are coded as 1s and 0s, repectively. I'll keep that for now, but will change in a bit
			genotypes_as_freqs = (variant.genotypes[::2] + variant.genotypes[1::2])/2
			
			freqs_per_pop = ([i.mean()  for i in chunks(genotypes_as_freqs,int(diploids_per_pop))])
			
			if sum(freqs_per_pop)/len(freqs_per_pop) < 0.1: continue
			
			freq_dict[variant.position] = freqs_per_pop

		for pair in itertools.product(freq_dict.keys(), repeat=2):
			pw_distances_ij = abs( pair[0] - pair[1] )
			if pw_distances_ij == 0: continue
			LD_output_between.write(",".join([str(int(pw_distances_ij)), str(scipy.stats.pearsonr(freq_dict[pair[0]], freq_dict[pair[1]])[0]**2) ] ) + "\n")

	LD_output_between.close()

	return



def summary_stats(tree_sequence_file, outputFileName):
	ts = pyslim.load(tree_sequence_file)
	getLD(ts, outputFileName)
	return

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--trees", 
			required = True,
			dest = "trees",
			type = str, 
			help = "Directory with the coalescent trees for the simulation")
	parser.add_argument("--number", 
			required = True,
			dest = "num",
			type = int, 
			help = "What number of LD output files do you want?")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "Give a suffix for the output files [ the generation and replicate # will be added")
	parser.add_argument("--gen", 
			required = False,
			dest = "gen",
			type = int, 
			help = "Give the generation number that you will analyse")
	args = parser.parse_args()
	
	count = 0
	for i in glob.glob(args.trees + "/*trees"):
		if count >= args.num:
			print("we have done it!")
			break
		name = i.split("/")[-1].split(".")
		print("analysing " + i)
		direc = "/".join(i.split("/")[0:-1])
		rep = name[0]

		if rep != "4_0": continue
		gen = name[1]
		if args.gen:
			if int(gen) != args.gen:
				continue
		else:
			pass
			
		outputFileName = rep + "." + gen + ".LD.csv"

		summStats = summary_stats(i, outputFileName)
		count += 1

main()
