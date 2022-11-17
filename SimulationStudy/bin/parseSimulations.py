import msprime, pyslim, random, allel, math
import pandas as pd
import scipy.stats
import numpy as np
import argparse, tskit
import matplotlib.pyplot as pyplot


## This function grabs a subSample of individuals from a subset of pops.
def getSample(tree_sequence, nPops, nInds):
	pops = {}
#	print("!!!")
	for i in pyslim.individuals_alive_at(tree_sequence, 0):
#		print(i)
#		print(tree_sequence.individual(i).metadata["subpopulation"])
		if tree_sequence.individual(i).metadata["subpopulation"] == 999: continue
		try:
			pops[tree_sequence.individual(i).metadata["subpopulation"]].append( i )
		except KeyError:
			pops[tree_sequence.individual(i).metadata["subpopulation"]] = [i]
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

def getPVE_stabilising(treeSequence, enviDictionary, Genome):
## Get the total haploid population size
	N = treeSequence.get_sample_size()

## Get the number of subpopulations - this is a pretty inelegant way of doing it...
	subPops =  set( [d for d in treeSequence.individual_populations if d != 999]  )

	environments = [enviDictionary[i] for i in sorted(subPops)]

#	individual_populations= [d for d in treeSequence.individual_populations if d != 999][:int(N/2)]

	windows = [[660000, 669999],
	[1320000, 1329999],
	[1980000, 1989999],
	[2640000, 2649999],
	[3300000, 3309999],
	[3960000, 3969999],
	[4620000, 4629999],
	[5280000, 5289999],
	[5940000, 5949999],
	[6600000, 6609999],
	[7260000, 7269999],
	[7920000, 7929999]]

# Define the windows ( if they weren't passed to the function directly )
# Add in windows for the three cartoon chromosomes at the end of the chromosome

##	Set up empty vectors for the metadata you want from the treeSequence (selection coefficients, positions etc.)
	selCoeffs = [ ]
	positionVectorCutter = [ ]
## At each of the variable in the tree, get the positions and the selection coefficients of each allele
## (remember that we haven't added neutral mutations to this tree yet!)
	for i in treeSequence.variants():
		positionVectorCutter.append(i.position)
		selCoeffVec = [0]
		for m in i.site.mutations:
			if len( m.metadata["mutation_list"] ) > 1:
				print("site has more than one mutation")
			selCoeffVec.append( m.metadata["mutation_list"][0]["selection_coeff"] )
		selCoeffs.append(np.array(selCoeffVec))

# If there are no SNPs, just return an np.empty()
	if len(selCoeffs) == 0:
		return np.empty(len(windows)), np.empty(len(windows))

# Make the (potentially) ragged array of SNP effects for each allele at each position
	selCoeffs_all = np.array(selCoeffs)

## This array is called the cutter because it is used to slice up the Window SNP stuff later.
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
#	C_p = np.var(C_p_vec) # Between population phenotypic variance
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
#	C_pg = C_pg_vec/C_p
	Cov_pg = Cov_pg_vec/Cov_p

	adapGenes = {}
	for w in range(len(windows)):
		print( Genome[Genome[0] == windows[w][0]+1].names.values[0] )
		focal_gene = Genome[Genome[0] == windows[w][0]+1].names.values[0]
		adapGenes[ focal_gene ] = Cov_pg[w]

		if Cov_pg[w] != 0:
# Mark the 3 genes up and downstream from the targets of selection
			for k in [1,2,3]:
				adapGenes[ "gene" + str(int(focal_gene[4:]) + k) ] = -99
				adapGenes[ "gene" + str(int(focal_gene[4:]) - k) ] = -99


	return adapGenes

def getAdapGenes(PVE, Genome):
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
		## All sites are biallelic (i.e. there are only two possible states, so if there is a mutation present whether it's a single allele or not, they:
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

#	for d in treeSequence.individuals():
#		print(d.metadata["subpopulation"])

#	individual_populations= [d for d in treeSequence.populations() if d != 999][:int(N/2)]
	individual_populations= [d.metadata["subpopulation"] for d in treeSequence.individuals()]

	envArray = np.array( [enviPops[d] for d in individual_populations] )

## At each of the variable in the tree, get the positions and the selection coefficients of each allele
## (remember that we haven't added neutral mutations to this tree yet!)

	geneDict = {}

	for i in treeSequence.variants():
		positionVectorCutter.append(i.position)
#		print(i.genotypes)

#		i.genotypes[i.genotypes == 2] = 1
		if sum(i.genotypes == 2) != 0:
			print("! multiallelic site")
#		## We do not bother with multiallelic sites
			continue

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


## A little function to write VCF lines
def VCF_line(ts_variant, previous_position):
	current_pos = round(ts_variant.position)
	if current_pos == round(previous_position):
		current_pos +=1
#	genotype_list_raw = [g for g in ts_variant.genotypes]

	genotype_list_raw = [k if k>=0 else 0 for k in ts_variant.genotypes]
	genotype_list = [str(genotype_list_raw[g])+"|"+str(genotype_list_raw[g+1]) for g in range(len(genotype_list_raw))[::2]]
	strt = "1	"+str(current_pos)+"	.	0	1	.	PASS	.	GT	"
	return(strt+"\t".join(genotype_list)+"\n", current_pos)

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
			type = float,
			help = "If these simulations are modelling directional selection, what was s?",
			default = 0.003)
	parser.add_argument("--environments",
			required = False,
			dest = "environments",
			type = str,
			help = "The file with environments in it - if not given, optima will be used")
	parser.add_argument("--demo",
			required = False,
			dest = "demo",
			action = "store_true",
			help = "Do you want to make a file with allele frequency info for demonstration purposes?")
	parser.add_argument("--VCF",
			required = False,
			dest = "VCF",
			action = "store_true",
			help = "Do you want to make a VCF file?")
	parser.add_argument("--mutMap",
			required = False,
			dest = "mutMap",
			action = "store_true",
			help = "Do you want to simulate varying mutation rates?")
	parser.add_argument("--mut_mod",
			required = False,
			dest = "mut_mod",
			type = float,
			help = "The factor to multiply the mutation rate by",
			default = 1.)
	args = parser.parse_args()

# Sample mutation rates as log10(mutation rates)
	if args.mutMap:
		mutRates =  10**(np.random.uniform(low = -9, high = -7.43, size = 1000) )
		boundaries = [ep for ep in range(0,10010000, 10000)]
		boundaries[-1] +=1
		mutRateMap = msprime.RateMap( position = boundaries, rate = mutRates)

	print("analysing",args.trees)

## Read in the tree sequence file
	if args.nPops == 0:
		if args.VCF:
			print("can only make VCFs for samples")
			return
		ts = pyslim.update( tskit.load(args.trees) )

	else:
		ts_raw = pyslim.update( tskit.load(args.trees) )
		ts = getSample(ts_raw, args.nPops, args.nInds)

## The simulations have a common genome structure,
## mimicking a group of species with syntenic genomes,
	genome = genomeMaker()
	genome['names'] = np.array(['gene'+str(i) for i in range(genome.shape[0])])

## Get the sample size
	N = ts.get_sample_size()

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
	for i in pyslim.individuals_alive_at(ts, 0):
		count += 1

		if ts.individual(i).metadata["subpopulation"] == 999: continue
		try:
			subPopDict[ts.individual(i).metadata["subpopulation"]].append( i )
		except KeyError:
			subPopDict[ts.individual(i).metadata["subpopulation"]] = [i]

	print(count,"individuals")
	envi_pops = [enviDict[i] for i in sorted(subPopDict.keys()) ]
	print(envi_pops)

	nPops = len( list(subPopDict.keys()) )
	diploids_per_pop = len( subPopDict[list(subPopDict.keys())[0]] )

	print( diploids_per_pop ,"diploid individuals per population")
	print( nPops ,"populations")


#	print(envi_pops, len(envi_pops))
	if args.bayPass:
		bayPassEnvs = open(args.output + ".bayPass.pc1", "w")
		for e in envi_pops:
			bayPassEnvs.write(str(e) + " ")
		bayPassEnvs.close()

## Figure out from the whole population, which genes are involved in local adaptation or not
	if not args.directional:
		if args.nPops == 0:
			adapGenes = getPVE_stabilising(ts , optimaDict, genome)
		else:
			adapGenes = getPVE_stabilising(ts_raw , optimaDict, genome)
	else:
		if args.nPops == 0:
			adapGenes = getPVE_Directional(ts , diploids_per_pop, optimaDict, genome)
		else:
			adapGenes = getPVE_Directional(ts_raw , diploids_per_pop, optimaDict, genome)

#	print(adapGenes)
#	for i in ts.individuals_alive_at(0):
#		print(ts.individual(i).id, ts.individual(i).population, enviDict[ts.individual(i).population])
## Uncomment for the time bing
#	if len( enviDict.keys() ) != len( subPopDict.keys() ):
#		print( "Uneven numbers of subpopulations in the simulation and the optima file that was given" )
#		return




## Set the mutation rate for neutral mutations
	mut_rate = 1e-8 * args.mut_mod ## HARD CODED
#	mut_rate = 1e-8/1000 ## HARD CODED

## If you specify an interval, the following will cut down the tree to focus sub-windows within each 10,000bp window
## This is a way of getting low recombination results on the cheap
## ...if you are into the whole brevity thing
	if args.intervals != 0:
		intervalDiff = (10000 - args.intervals)/2
		if args.intervals == 1:
			print("I have not put in support for single base pair intervals")
			return
		genome["interval_start"] = genome[0] + intervalDiff
		genome["interval_end"] = genome[1] - intervalDiff
		intervalArray = np.array(genome[["interval_start","interval_end"]])
		print( genome )
		ts = ts.keep_intervals( intervalArray )

	## Reset the mutation rate for neutral mutations to acheive the same net mutation rate
		mut_rate = 1e-8 * (10000/args.intervals)
 # HARD CODED, SO CHANGE IF NECESSARY

#	print("building VCF")
#
	if args.VCF:
 		## Writing VCFs requires a fully coalesced population
		if args.mutMap:
			sprinkled = msprime.mutate(ts, rate= mutRateMap, keep=True, random_seed = 12345)
		else:
			sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True, random_seed = 12345)
#
		ind_names = []
#
		envi_output = open(args.output+"_individualEnvs.txt", "w")
#
		for ev in envi_pops:
				for pp in range(args.nInds):
					envi_output.write(str(ev) + "\n")
		envi_output.close()
		for POP in range(args.nPops):
			for IND in range(args.nInds):
				ind_names.append("p"+str(POP)+"_"+"i"+str(IND)+"_")
#
#
		vcf_file = open(args.output + ".vcf", "w")
		vcf_file.write("##fileformat=VCFv4.2\n")
		vcf_file.write("##source=tskit 0.4.1\n")
		vcf_file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
		vcf_file.write('##contig=<ID=1,length=10000001>\n')
		vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">")\n')
		vcf_file.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	'+"\t".join(ind_names)+"\n")
		prev_pos = 0
		for variant in sprinkled.variants():
#
 #			print(variant)
			output_line, current_position = VCF_line(variant, prev_pos)
			vcf_file.write(output_line)
# #			print(variant.genotypes)
# #			print(VCF_line(variant, prev_pos))
# #			print()
			prev_pos = current_position
# #			break
		vcf_file.close()
#		return

## Sprinkle neutral mutations onto the coalescent tree
	print('Sprinkling mutations onto trees')

	if args.mutMap:
		sprinkled = msprime.mutate(ts, rate= mutRateMap, keep=True, random_seed = 12345)
	else:
		sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True, random_seed = 12345)

	# if args.VCF:
	#
	# 	ind_names = []
	#
	# 	for POP in range(args.nPops):
	# 		for IND in range(args.nInds):
	# 			ind_names.append("p"+str(POP)+"_"+"i"+str(IND)+"_")
	#
	# 	print("writing VCF")
	# 	with open(args.output+".vcf", "w") as vcf_file:
	# 		sprinkled.write_vcf(vcf_file, position_transform = "legacy")

	print('extracting segregating sites from trees...')
# Make a little dict that will be populated by the genes that contribute to LA

	selSites = 0
	count = 0

	outputDF = open( args.output +'.csv' ,'w' )

	if args.demo:
		print("I'm not performing GEA analysis, I'm make an allele frequency SNP table for demonstration purposes")
		envOutput = open("environments."+args.output+".csv", "w")
		envOutput.write( ",".join( map(str, envi_pops) ) + "\n" )
		envOutput.close()
	else:
		print( "Not demo!" )
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
		gene = findGene(variant.position, genome)

		if not gene:
			continue

		elif gene == -99:
			print(gene)
			print('FAILED')
			return
	#	print(variant.position)

#		md = pyslim.decode_mutation(variant.site.mutations[0].metadata)
		if args.directional:
			if len(variant.site.mutations[0].metadata["mutation_list"] ) > 0:
			#print( len(variant.site.mutations[0].metadata["mutation_list"]) )
				selCoeff = args.s
			else:
				selCoeff = 0
		else:
			variant.site.mutations[0].metadata["mutation_list"]

#			md = pyslim.default_slim_metadata(variant.site.mutations[0])
			if len(variant.site.mutations[0].metadata["mutation_list"] ) > 0:
				selCoeff =  variant.site.mutations[0].metadata["mutation_list"][0]["selection_coeff"]
#				selCoeff = md[0].selection_coeff
			else:
				selCoeff = 0
#		print()

#		continue

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

		if args.s != 0:

			try:
				LA = adapGenes[gene[0]]
			except KeyError:
				LA = 0
		else:
			LA = 0

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

		if args.demo == True:
#			print("DEMO")
			outline = [ "chr" + str(1 + math.floor(variant.position/2000000)),
					int(variant.position),
					gene[0] ] +  freqs_per_pop
#			print(outline)
			outputDF.write(",".join(map(str, outline)) + "\n")

			continue



#		geno_spearman = scipy.stats.spearmanr(envi, genotypes_as_freqs)
#		geno_k_tau, geno_k_tau_p_value = scipy.stats.kendalltau(envi, genotypes_as_freqs)
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
					pbar,
					pbar_qbar,
					maf,
					LA]

		outputDF.write(",".join(map(str, outline)) + "\n")

	outputDF.close()

	if args.bayPass:
		bayPassConfig.close()
main()
