## A script to extract sub-trees for directionally selected sites to get a better understanding of the effect of pbar*q_bar weighting

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
		sampled_pops = [153, 173, 143,  18,  53,  60, 152, 148,  15,  47, 166,  51, 121,  50, 149, 101, 111,  76,   1, 110]
	elif nPops == 10:
## Fix the sample so that the same pops are sampled across replicates
		sampled_pops = [1,  51, 149, 110, 148,  15, 101, 153,  47, 166]
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




def getPVE_Directional(treeSequence, diploidsPerPop, enviPops, genome, s = 0.005):
## Get the total haploid population size
	N = treeSequence.get_sample_size()
#	print(N)
## Get the number of subpopulations - this is a pretty inelegant way of doing it... 
#	subPops =  set( [d for d in treeSequence.individual_populations if d != 999]  )

##	Set up empty vectors for the metadata you want from the treeSequence (selection coefficients, positions etc.)

	selCoeffs = [ ] 

	positionVectorCutter = [ ] 

	genotypeFrequencies = [ ]

#	print( [d for d in treeSequence.individual_populations], len([d for d in treeSequence.individual_populations]) )
#	print( [enviPops[d-1] for d in treeSequence.individual_populations if d != 999] )

	individual_populations= [d for d in treeSequence.individual_populations if d != 999][:int(N/2)]

	envArray = np.array( [enviPops[d] for d in individual_populations] )

## At each of the variable in the tree, get the positions and the selection coefficients of each allele
## (remember that we haven't added neutral mutations to this tree yet!)

	geneDict = {}

	for i in treeSequence.variants():
		positionVectorCutter.append(i.position)
		i.genotypes[i.genotypes == 2] = 1

		if i.genotypes.sum()/len( i.genotypes ) >= 1:
		## We do not bother with fixed mutations in the directional selection model
			continue

		genotypes_as_freqs = (i.genotypes[::2] + i.genotypes[1::2])/2

		print( genotypes_as_freqs )

		relevantGene =  findGene(i.position, genome) 
		if len(relevantGene) > 1:
			print("WHAT! How can a SNP be present within 2 or more genes?")

## Calculate fitness contribution of each variant to population fitness 
		fitnessContribution = genotypes_as_freqs*envArray*s
## Calculate the covariance of fitness contributions and environments
		fitnessEnvCovariance = np.cov( fitnessContribution, envArray )[0,1]
		geneDict[ relevantGene[0] ] = fitnessEnvCovariance
		for k in [1,2,3]:
			geneDict[ "gene" + str(int( relevantGene[0][4:] ) + k) ] = -99
			geneDict[ "gene" + str(int( relevantGene[0][4:] ) - k) ] = -99


	return geneDict
	
	
	
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
	parser.add_argument("--nPops", 
			required = False,
			dest = "nPops",
			type = int, 
			help = "The number of subPops to downsample to",
			default = 0)
	parser.add_argument("--nInds", 
			required = False,
			dest = "nInds",
			type = int, 
			help = "The number of individuals to grab from each subPop",
			default = 0)
	parser.add_argument("--bayPass", 
			required = False,
			dest = "bayPass",
			action = "store_true",
			help = "Do you want to spit out a BayPass config?")
	parser.add_argument("--directional", 
			required = False,
			dest = "directional",
			action = "store_true",
			help = "Are these simulations modelling directional selection?")
	parser.add_argument("--intervals", 
			required = False,
			dest = "intervals",
			type = int, 
			help = "The width of the intervals you want to analyse",
			default = 0)
	parser.add_argument("--s", 
			required = False,
			dest = "s",
			type = int,
			help = "Are these simulations modelling directional selection?")
	parser.add_argument("--environments", 
			required = False,
			dest = "environments",
			type = str, 
			help = "The file with environments in it - if not given, optima will be used")
	args = parser.parse_args()

	print("analysing",args.trees)
	
## Read in the tree sequence file
	if args.nPops == 0:
		ts = pyslim.load(args.trees)
	else:
		ts_raw = pyslim.load(args.trees)
		ts = getSample(ts_raw, args.nPops, args.nInds)


## Make a dict of the environments from the optima and (optional) enviroment files
	if args.environments == None:
		enviDict = {}
		count = 0
		with open(args.optima) as file:
			for line in file:
				enviDict[ count ] = float( line.strip() )
				count += 1
		optimaDict  = enviDict.copy()
	else:
		enviDict = {}
		count = 0
		with open(args.environments) as file:
			for line in file:
				enviDict[ count ] = float( line.strip() )
				count += 1
		optimaDict = {}
		count = 0
		with open(args.optima) as file:
			for line in file:
				optimaDict[ count ] = float( line.strip() )
				count += 1
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
	print(envi_pops)
	


	to_keep = []
	for v in ts.variants():
		if v.genotypes.sum()/(args.nPops*args.nInds) < 0.3: continue
		to_keep.append( [ v.position-0.01, v.position + 0.01] )
		
	ts = ts.keep_intervals( np.array(to_keep) )
	
#	for t in ts_edit.trees():
#		print(t.draw(format = "unicode"))

	if len( [q.position for q in ts.variants()]  )*2+1 != len( [q.interval for q in ts.trees()]  ):
		print("WTF, this is error 1")

	
## Get the sample size
	N = ts.get_sample_size()



	nPops = len( list(subPopDict.keys()) )
	diploids_per_pop = len( subPopDict[list(subPopDict.keys())[0]] )

	print( diploids_per_pop ,"diploid individuals per population")
	print( nPops ,"populations")


	print(envi_pops, len(envi_pops))


## Set the mutation rate for neutral mutations 
	mut_rate = 5 ## HARD CODED

## Sprinkle neutral mutations onto the coalescent tree
	print('Sprinkling mutations onto trees')

	sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True, random_seed = 12345)

	print('extracting segregating sites from trees...')
# Make a little dict that will be populated by the genes that contribute to LA

	selSites = 0
	count = 0

	outputDF = open( args.output +'.csv' ,'w' )
	header = [ "position",
				"gene",
				"selCoeff", 
				"geno_spearman_correlation",
				"geno_spearman_correlation_2",
				"geno_spearman_pvalue",
				"geno_k_tau",
				"geno_k_tau_p_value",
				"pbar",
				"pbar_qbar",
				"maf",
				"LA" ]
	outputDF.write(",".join(header) + "\n")
	
	if args.bayPass:
		
		bayPassConfig = open( args.output + ".bayPass.txt", "w")
	
## Iterating over all variants in the population we now perform the actual GEA	

	for variant in sprinkled.variants():
#		if variant.position > 1e6:
#			continue
		if variant.num_alleles > 2:
			continue

		md = pyslim.decode_mutation(variant.site.mutations[0].metadata)

		if len(md) > 0:
			selCoeff = md[0].selection_coeff
		else:
			selCoeff = 0

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
		
		if variant.genotypes.sum() == N: continue # Ignore fixed mutations

		freqs_per_pop = ([i.mean()  for i in chunks(genotypes_as_freqs,int(diploids_per_pop))])


		pbar = sum(freqs_per_pop)/len( freqs_per_pop )

		LA = -99
		
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

		geno_spearman = scipy.stats.spearmanr(envi_pops, freqs_per_pop)
		geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(envi_pops, freqs_per_pop)
		
		pop_spearman = scipy.stats.spearmanr(envi_pops, freqs_per_pop)
		pop_k_tau, pop_k_tau_p_value = scipy.stats.kendalltau(envi_pops, freqs_per_pop)
		
		if variant.position == int(variant.position):
			selCoeff = 0.003
		else:
			selCoeff = 0.0
		gene = int(round(variant.position/10)*10)
			
		outline = [ variant.position, 
					gene,
					selCoeff, 
					geno_spearman.correlation,
					geno_spearman.correlation**2,
					geno_spearman.pvalue,
					geno_k_tau,
					geno_k_tau_p_value,
					pbar,
					pbar_qbar,
					maf,
					LA]
					
		outputDF.write(",".join(map(str, outline)) + "\n")

	outputDF.close()

	if args.bayPass:
		bayPassConfig.close()
main()
