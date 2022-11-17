import msprime, pyslim, random, allel
import pandas as pd
import scipy.stats
import numpy as np
import argparse, tskit
import matplotlib.pyplot as pyplot


# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
# For item i in a range that is a list of length l,
	for i in range(0, len(l), n):
		# Create an index range for l of n items:
		yield l[i:i+n]


## This function grabs a subSample of individuals from a subset of pops.
def getSample(tree_sequence, nPops, nInds):
	pops = {} 
	for i in tree_sequence.individuals_alive_at(0):
		if tree_sequence.individual(i).population == 999: continue

		try: 
			pops[tree_sequence.individual(i).population].append( i )
		except KeyError:
			pops[tree_sequence.individual(i).population] = [i]
#	print([i for i in pops.keys()])
# Let's take a sample of 40 pops and 20 diploids from each and make a new tree from them
	sampled_indivs = []

	if nPops == 40:
## Fix the sample so that the same pops are sampled across replicates
		sampled_pops = [119, 15, 148, 106, 104, 111, 18, 76, 50, 173, 166, 177, 101, 1, 153, 35, 149, 137, 40, 121, 110, 124, 89, 113, 98, 115, 51, 47, 96, 8, 120, 152, 91, 60, 31, 117, 143, 53, 49, 81]
## Fix the sample so that the same pops are sampled across replicates
	elif nPops == 20:
		sampled_pops = [147,   3,  59,  27,  58, 122, 143, 160,  75,  26,  20,  94,  34, 165, 161,  90,   2,  92, 151,  64]
	elif nPops == 10:
## Fix the sample so that the same pops are sampled across replicates
		sampled_pops = [91,  99, 111, 180, 166,  15,  46, 146, 154,  98]
	else:
		## The following would randomly sample pops
		sampled_pops = np.random.choice( list( pops.keys() ), nPops, replace = False)

	print(sampled_pops, len(sampled_pops))
#	for i in np.random.choice( list( pops.keys() ), nPops, replace = False):
	for i in sampled_pops:
		for j in np.random.choice( np.array(pops[i]), nInds , replace = False):
			sampled_indivs.extend( tree_sequence.individual(j).nodes )
	sampledTree = tree_sequence.simplify(sorted(sampled_indivs))
	return sampledTree


def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--trees", 
			required = True,
			dest = "trees",
			type = str, 
			help = "Coalescent trees for the simulation")
	parser.add_argument("--optima", 
			required = True,
			dest = "optima",
			type = str, 
			help = "The file of phenotypic optima for the simulation")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "What name do you want to give to the output file? [DON'T GIVE FILE EXTENSION, THE SCRIPT MAKES TWO OUTPUT FILES")
	args = parser.parse_args()

	print("analysing",args.trees)
	
## Read in the tree sequence file
	ts_raw = pyslim.load(args.trees)
	nPops = 20
	nInds = 100
	
	ts = getSample(ts_raw, nPops, nInds)

## The simulations have a common genome structure, 
## mimicking a group of species with syntenic genomes,

## Get the sample size
	N = ts.get_sample_size()

## Make a dict of the environments from the optima and (optional) enviroment files
	enviDict = {}
	count = 0
	with open(args.optima) as file:
		for line in file:
			enviDict[ count ] = float( line.strip() )
			count += 1
	optimaDict  = enviDict.copy()


## Let's make a dict of all the demes and the individuals from the tree being analysed
	subPopDict = {} 
	count = 0
	for i in ts.individuals_alive_at(0):
		count += 1
		if ts.individual(i).population == 999: continue
		try: 
			subPopDict[ts.individual(i).population].append( i )
		except KeyError:
			subPopDict[ts.individual(i).population] = [i]

	print(count,"individuals")
	envi_pops = [enviDict[i] for i in sorted(subPopDict.keys()) ]
	
	nPops = len( list(subPopDict.keys()) )
	diploids_per_pop = len( subPopDict[list(subPopDict.keys())[0]] )

	print( diploids_per_pop ,"diploid individuals per population")
	print( nPops ,"populations")


	print(envi_pops, len(envi_pops))


## Set the mutation rate for neutral mutations 
	mut_rate = 1e-8 ## HARD CODED

## Sprinkle neutral mutations onto the coalescent tree - specify the random seed to make the analysis reproducible
	print('Sprinkling mutations onto trees')

	sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True, random_seed = 12345)

	print('extracting segregating sites from trees...')

# Make a little dict that will be populated by the genes that contribute to LA

	output = open(args.output, "w")
## Iterating over all variants in the population we now perform the actual GEA	

	for variant in sprinkled.variants():
		if variant.position < 8000000:continue
		print(variant.position)
		if variant.num_alleles > 2:
			continue


		count +=1
		if count % 1000 == 0:
			print('extracted',count,'neutral sites')

## Calculating the genotype freqs this way assumes that the ref and alt are coded as 1s and 0s, repectively. I'll keep that for now, but will change in a bit

		if 2 in variant.genotypes:
			print(variant)
			print(list(variant.genotypes))
			print("!!!!!!!!!!!!!!")
			close()

		genotypes_as_freqs = (variant.genotypes[::2] + variant.genotypes[1::2])/2
		
		print(genotypes_as_freqs.shape)
		
		if variant.genotypes.sum() == N: continue # Ignore fixed mutations

		freqs_per_pop = ([i.mean()  for i in chunks(genotypes_as_freqs,int(diploids_per_pop))])


		pbar = sum(freqs_per_pop)/len( freqs_per_pop )


		if pbar == 1:continue ## This can be triggered by mutation stacking and infinite sites funniness, fix this part of the script
		if pbar > 1:
			print(variant)
			print("WTF, what is going on here: pbar > 1")
			print(genotypes)
			continue

		pbar_qbar = pbar * ( 1 - pbar )
		maf = min(pbar, 1-pbar)
		if pbar_qbar <0:
			print("WTF, what is going on here: pbar_qbar < 1")
			print(genotypes)
			continue

		
		true_cor_geno_spearman = scipy.stats.spearmanr(envi_pops, freqs_per_pop)
		true_cor_geno_k_tau, true_cor_geno_k_tau_p_value = scipy.stats.kendalltau(envi_pops, freqs_per_pop)


		sample_reps = 100
		correlation_estimates = []
		for q in range( sample_reps ):
			freqs_per_pop_sample = [ np.random.choice(i, 30).mean() for i in chunks(genotypes_as_freqs,int(diploids_per_pop)) ]
			if sum(freqs_per_pop_sample) == 0 or sum(freqs_per_pop_sample) == 1: continue
			sample_cor_geno_spearman = scipy.stats.spearmanr(envi_pops, freqs_per_pop_sample)
			correlation_estimates.append(sample_cor_geno_spearman.correlation )
		correlation_estimates = np.array( correlation_estimates )
		output.write( ",".join([str(np.var(correlation_estimates)), str(maf), str(pbar_qbar)] ) + "\n")


	output.close()

main()



