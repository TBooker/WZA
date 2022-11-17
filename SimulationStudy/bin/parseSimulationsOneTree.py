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
#	sampled_pops = np.random.choice( list( pops.keys() ), nPops, replace = False)
	sampled_pops = [119, 15, 148, 106, 104, 111, 18, 76, 50, 173, 166, 177, 101, 1, 153, 35, 149, 137, 40, 121, 110, 124, 89, 113, 98, 115, 51, 47, 96, 8, 120, 152, 91, 60, 31, 117, 143, 53, 49, 81]

	print(sampled_pops)
#	for i in np.random.choice( list( pops.keys() ), nPops, replace = False):
	for i in sampled_pops:
		for j in np.random.choice( np.array(pops[i]), nInds , replace = False):
			sampled_indivs.extend( tree_sequence.individual(j).nodes )
	sampledTree = tree_sequence.simplify(sorted(sampled_indivs))
	return sampledTree

def getPVE(treeSequence, environments, windows = None):
## Get the total haploid population size
	N = treeSequence.get_sample_size()
	
## Get the number of subpopulations - this is a pretty inelegant way of doing it... 
	subPops =  set( [d for d in treeSequence.individual_populations if d != 999]  )

# Define the windows ( if they weren't passed to the function directly )
# Add in windows for the three cartoon chromosomes at the end of the chromosome

	if windows == None:
		windows = [[49500+ (i*100000),59499 + (i*100000)] for i in range(10)]
	windows += [[1e6,1e6+999],[1e6+1000,1e6+1999],[1e6+2000,1e6+2999]]

#	if windows == None:
#		windows = [[49500+ (i*100000),59500 + (i*100000)] for i in range(20)]

##	Set up empty vectors for the metadata you want from the treeSequence (selection coefficients, positions etc.)
	selCoeffs = [ ] 
	positionVectorCutter = [ ] 
## At each of the variable in the tree, get the positions and the selection coefficients of each allele
## (remember that we haven't added neutral mutations to this tree yet!)
	for i in treeSequence.variants():
		positionVectorCutter.append(i.position)
		selCoeffVec = [0]
		for m in i.site.mutations:
			if len( m.metadata ) > 1:
				print("site has more than one mutation")
			selCoeffVec.append( m.metadata[0].selection_coeff )
		selCoeffs.append(np.array(selCoeffVec))
		
# If there are no SNPs, just return an np.empty()
	if len(selCoeffs) == 0:
		return np.empty(len(windows)), np.empty(len(windows))
		
# Make the (potentially) ragged array of SNP effects for each allele at each position
	selCoeffs_all = np.array(selCoeffs)

## This array is called the cutter because it is used to slice up the Window SNP stuff later. Plus The Cutter is a cool song by E&TBM
	positionVectorCutter = np.array(positionVectorCutter)

## Defining numpy arrays as I have done here is probably looked down upon by some, 
## however, getting these arrays is by no means the limiting step, so meh

## Grab the genotype matrix for the SNPs
	msprime_genotype_matrix = treeSequence.genotype_matrix()

# Convert msprime's haplotype matrix into a sci-kit allel genotype matrix by merging chromosomes
	haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )

# Convert from a haploid to a genotype matrix (with a length 2 array for each element) 
	genotype_array = np.array(haplotype_array.to_genotypes(ploidy=2))

# Initialise an empty array that will hold the phenotypic contribution of each variable site to an individual's phenotype
	alpha_array = np.empty( shape = genotype_array.shape )

	for k in range( genotype_array.shape[0]):
		alpha_array[k] = selCoeffs_all[k][ genotype_array[k] ]

# Now sum the alpha for the two genotypes
	alpha = alpha_array.sum(axis = 2 ).transpose()
# Get the population size per deme (the // operator is integer division in Python3.x)
	sizePerDeme = alpha.shape[0]//len(subPops)

## The code below can be used to get the mean for each population
#	count = 0
#	for i in alpha.sum(axis = 1 ):
#		print(i)
#		count+=1
#		if count ==20:break
#	phens = alpha.sum(axis = 1 )
#	means = []
#
#	for x in range( len (subPops) ):
#		temp = phens[x*sizePerDeme:(x+1)*sizePerDeme,]
#		means.append(sum(temp)/len(temp))
#	print(means)
#	sys.exit()

# Sanity check: the number of sites in the alpha array (i.e. the second dimension) should
# be equal to the number of selected sites
	if alpha.shape[1] != len( selCoeffs ):
		print("there is a mismatch between the size of the genotype matrix and the SNP effect vector")
		print(genotype_array.shape[1], len(selCoeffs) )
		sys.exit()
		
# Initialise an empty vector to store the C_p data in (The vector of mean phenotypes for each population)
# We take the genotype matrix for all individuals and sum across sites to get individual phenotypes
# The average of these across populations is then used to get the phenotypic variance 
	C_p_vec = np.empty(len(subPops))
	
	for i in range( len (subPops) ):
		subMatrix =  alpha[i*sizePerDeme:(i+1)*sizePerDeme,:]
		C_p_vec[i] = np.mean(subMatrix.sum(axis = 1))
#		print( subMatrix.sum(axis = 1)) 
# This gets the total additive genetic variance for all SNPs
	C_p = np.var(C_p_vec) # Between population phenotypic variance
	Cov_p = np.cov(C_p_vec, environments)[0,1] # Between population  covariance of phenotype and environment

#	print("here's the population-level variance", C_p)
#	print("here's the population-level covariance", np.cov(C_p_vec, environments))

# Initialise an empty vector to store the C_pg data in for each gene
	C_pg_vec = np.empty(len(windows))
	Cov_pg_vec = np.empty(len(windows))

	for win in range(len(windows)):
		# Get the SNPs in the window
		relevantSNPs = (positionVectorCutter >=windows[win][0])&(positionVectorCutter <=windows[win][1])
		# Initialise an empty array
		C_pdg = np.empty(len(subPops))
		
		for i in range( len (subPops) ):
			subMatrix =  alpha[i*sizePerDeme:(i+1)*sizePerDeme,relevantSNPs]
#			if len(selCoeffs[relevantSNPs]) > 0:
#				print("\n"+str(i))
#				print(selCoeffs[relevantSNPs])
#			print(subMatrix, subMatrix.shape)
			C_pdg[i] =  np.mean( subMatrix.sum(axis = 1) )
#			print(  np.matmul( subMatrix, selCoeffs[relevantSNPs] ))
		C_pg_vec[win] = np.var(C_pdg)
		Cov_pg_vec[win] = np.cov(C_pdg, environments)[0,1]
#		print("here's the variance",np.var(C_pdg))
#		print("here's the covariance",np.cov(C_pdg, environments), np.cov(C_pdg, environments)[0,1])
	C_pg = C_pg_vec/C_p
	Cov_pg = Cov_pg_vec/Cov_p
	return C_pg, Cov_pg

def getAdapGenes(PVE):
	adapGenes = {}
	for i in range(0,10):
		adapGenes["gene" + str( (i *5) + 2 ) ] = PVE[i]
	return adapGenes

def getAdapGenesDirectional(treeSequence, diploidsPerPop, enviPops, genome):
	subPops =  set( [d for d in treeSequence.individual_populations if d != 999]  )
	C_p_vec = np.empty(len(subPops))
	
	gene_dict = {}
	for g in list(genome.names):
		gene_dict[g] = 1 ## All genes get a p-value of 1 unles they 

	position_vector_cutter = []
	correlation_vector = []
	for variant in treeSequence.variants():
#		if variant.num_alleles > 2:
#			continue
## Calculating the genotype freqs this way assumes that the ref and alt are coded as 1s and 0s, repectively. I'll keep that for now, but will change in a bit
		## A quick and dirt hack for testing purposes:
		variant.genotypes[variant.genotypes == 2] = 1

		if 2 in variant.genotypes:
			print(variant)
			print(list(variant.genotypes))
			print("!!!!!!!!!!!!!!")
			return
		relevantGene =  findGene(variant.position, genome) 
		if len(relevantGene) > 1:
			print("WHAT! How can a SNP be present within 2 or more genes?")
	
		genotypes_as_freqs = (variant.genotypes[::2] + variant.genotypes[1::2])/2
		
#		if variant.genotypes.sum() == N: continue # Ignore fixed mutations

		freqs_per_pop = ([i.mean()  for i in chunks(genotypes_as_freqs,int(diploidsPerPop))])
		pbar = sum(freqs_per_pop)/len( freqs_per_pop )
		if pbar == 1:continue
		maf = min(pbar, 1-pbar)

		geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(enviPops, freqs_per_pop)

#		position_vector_cutter.append(variant.position)

		gene_dict[relevantGene[0]] = geno_k_tau_p_value

#	position_vector_cutter = np.array(position_vector_cutter)
#	correlation_vector = np.array(correlation_vector)

#	correlations_per_gene = np.zeros(len(windows))

#	for win in range(len(windows)):
# Get the SNPs in each gene window
#		relevantSNPs = (position_vector_cutter >=windows[win][0])&(position_vector_cutter <=windows[win][1])
#		if relevantSNPs.sum()==0:continue
#		correlations_per_gene[win] = correlation_vector[relevantSNPs][0]
#	print(correlations_per_gene)
	return gene_dict


def genomeMaker_DEP():

#	geneStarts = np.array( [ i for i in range(0,1000000,10000) ][1::2] ) - 500
	geneStarts = np.array( [ i for i in range(0,1000000,5000) ][1::4] ) 
#	geneStarts = np.array( [ i for i in range(0,2000000,5000) ][1::4] ) 

#	geneEnds = geneStarts + 999
	geneEnds = geneStarts + 9999

	return( pd.DataFrame([geneStarts, geneEnds]).transpose() )


def genomeMaker():
	geneStarts = np.array( [ i for i in range(1,10000000,10000) ] ) 
	geneEnds = geneStarts + 9999
	return( pd.DataFrame([geneStarts, geneEnds]).transpose() )




# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
# For item i in a range that is a list of length l,
	for i in range(0, len(l), n):
		# Create an index range for l of n items:
		yield l[i:i+n]


## This function takes a list of gene intervals and a SNP position and for SNPs 
## inside the genes, returns the gene name
def findGene(pos, genes):
	geneID = (genes[0] <= pos) & (genes[1] >= pos) 
	myGene = list(genes[geneID].names.values)
	if len(myGene) > 1:
		print('Something got really messed up, investigate')
		return -99
	elif len(myGene) == 1:
		return myGene
	elif len(myGene) == 0:
		return None
	return

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="Perform GEA on SNPs present within a single Tree from each gene in the simulated genome")
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

	args = parser.parse_args()

	print("analysing",args.trees)
	
	if args.nPops == 0:
		ts = pyslim.load(args.trees)
	else:
		ts_raw = pyslim.load(args.trees)
		ts = getSample(ts_raw, args.nPops, args.nInds)



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


	if args.bayPass:
		bayPassEnvs = open(args.output + ".bayPass.pc1", "w")
		for e in envi_pops:
			bayPassEnvs.write(str(e) + " ")
		bayPassEnvs.close()

	nPops = len( list(subPopDict.keys()) )

	diploids_per_pop = len( subPopDict[list(subPopDict.keys())[0]] )

	print( diploids_per_pop ,"diploid individuals per population")
	print( nPops ,"populations")

## The simulations have a common genome structure, 
## mimicking a group of species with syntenic genomes,
	genome = genomeMaker()
	genome['names'] = np.array(['gene'+str(i) for i in range(genome.shape[0])])

	outputDF = open( args.output +'.csv' ,'w' )
	header = [ "position",
				"gene",
				"selCoeff", 
				"geno_spearman_correlation",
				"geno_spearman_correlation_2",
				"geno_spearman_pvalue",
				"geno_k_tau",
				"geno_k_tau_p_value",
				"pop_spearman_correlation",
				"pop_spearman_correlation_2",
				"pop_spearman_pvalue",
				"pop_k_tau",
				"pop_k_tau_p_value",
				"pbar",
				"pbar_q_bar",
				"maf",
				"LA" ]
	outputDF.write(",".join(header) + "\n")
		
	if not args.directional:
		winVars, winCovars = getPVE(ts , envi_pops)
		print(winVars)
		adapGenes = getAdapGenes(winCovars)
	else:
		adapGenes = getAdapGenesDirectional(ts , diploids_per_pop, envi_pops, genome)
		print(adapGenes)

	pos_ts = ts.keep_intervals(np.array([[pos,pos+1]]))


	print( np.array([[1,1]]) )
	print( np.array([[1,1]]).shape )
	
	for pos in np.array(genome[0]) + 4997:

#		pos_tree = tskit.TreeSequence( ts.at(pos) )

	## Set the mutation rate for neutral mutations 

		mut_rate = 1e-8  *10000 # HARD CODED, SO CHANGE IF NECESSARY

	## Sprinkle mutations onto the coalescent tree
		print('Sprinkling mutations onto trees')
		sprinkled = msprime.mutate(pos_tree, rate= mut_rate, keep=True)
		
	##	Calculating the PVE for each of the focal genes


	## Get the allele frequencies for all segregating sites at QTL
		print('extracting segregating sites from trees...')
		selSites = 0


		count = 0

		if args.bayPass:
			
			bayPassConfig = open( args.output + ".bayPass.txt", "w")
		
		for variant in sprinkled.variants():
	#		if variant.position > 1e6:
	#			continue
			if variant.num_alleles > 2:
				continue
			gene = findGene(variant.position, genome)

			if not gene:
				continue

			elif gene == -99:
				print(gene)
				print('FAILED')
				return
		#	print(variant.position)

			md = pyslim.decode_mutation(variant.site.mutations[0].metadata)
			if len(md) > 0:
				selCoeff = md[0].selection_coeff
			else:
				selCoeff = 0

			count +=1
	#		if count == 2:break
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

			if args.bayPass:
				counts_per_pop = ([int(i.sum()*2)  for i in chunks(genotypes_as_freqs,int(diploids_per_pop))])
				bayPassConfig.write(" ")
				for c in counts_per_pop:
					bayPassConfig.write(str(c) +" " + str(diploids_per_pop*2 - c)+" ")
				bayPassConfig.write("\n")
					
			
			try:
				LA = adapGenes[gene[0]]
			except KeyError:
				LA = 0	


			pbar = sum(freqs_per_pop)/len( freqs_per_pop )

			if pbar == 1:continue ## This can be triggered by mutation stacking and infinite sites funniness, fix this part of the script
			if pbar > 1:
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

			outline = [ variant.position, 
						gene[0], 
						selCoeff, 
						geno_spearman.correlation,
						geno_spearman.correlation**2,
						geno_spearman.pvalue,
						geno_k_tau,
						geno_k_tau_p_value,
						pop_spearman.correlation,
						pop_spearman.correlation**2,
						pop_spearman.pvalue,
						pop_k_tau,
						pop_k_tau_p_value,
						pbar,
						pbar_qbar,
						maf,
						LA]
			outputDF.write(",".join(map(str, outline)) + "\n")
		outputDF.close()
		if args.bayPass:
			bayPassConfig.close()

if __name__ == "__main__":
	main()
