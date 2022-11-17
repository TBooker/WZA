import msprime, pyslim, random, allel
import pandas as pd
import scipy.stats
import numpy as np
import argparse, tskit
import matplotlib.pyplot as pyplot



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
	for i in np.random.choice( list( pops.keys() ), nPops, replace = False):
		for j in np.random.choice( np.array(pops[i]), nInds , replace = False):
			sampled_indivs.extend( tree_sequence.individual(j).nodes )
	sampledTree = tree_sequence.simplify(sorted(sampled_indivs))
	return sampledTree

## This script just does some preliminary tests on the simulation data

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

	args = parser.parse_args()

	print("analysing",args.trees)
	
	ts = pyslim.load(args.trees)



	N = ts.get_sample_size()

	alive = ts.individuals_alive_at(0)

## Make a dict of the environments
	enviDict = {}
	count = 0
	with open(args.optima) as file:
		for line in file:
			enviDict[ count ] = int( line.strip() )
			count += 1

## Let's make a dict of all the demes and the individuals from each

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


## Set the mutation rate for neutral mutations 
	mut_rate = 1e-8 # HARD CODED, SO CHANGE IF NECESSARY

## Sprinkle mutations onto the coalescent tree
	print('Sprinkling mutations onto trees')
	sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True)
	

## Get the allele frequencies for all segregating sites at QTL
	print('extracting segregating sites from trees...')
# Make a little dict that will be populated by the genes that contribute to LA
#	adapGenes = {}	
	selSites = 0
## Iterate over all variants present on the tree from SLiM
	for vs in ts.variants():
		print(vs)
		print( vs.genotypes.sum()/len(vs.genotypes) )

#	for tree in ts.trees():
#		print(tree.interval, tree.interval[1] - tree.interval[0])

	sampled_trees = getSample(ts, 1, 10)
	
	count = 0
	for tree in sampled_trees.trees():
		if int(tree.interval[0])%10000 == 0 or int(tree.interval[1])%10000 == 0:
			count +=1
			print(tree.interval,int(tree.interval[0])%10000, int(tree.interval[0]), int(tree.interval[1]))
	print(count)

main()
