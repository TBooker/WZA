import msprime, pyslim, random
import pandas as pd
import scipy.stats
import numpy as np
import argparse, tskit
import matplotlib.pyplot as pyplot

def genomeMaker(genome_size = 20000000, gene_number = 4000):
	geneStarts = np.array( [i for i in range(1, genome_size, int(genome_size/gene_number))] )
	geneEnds = geneStarts + 99
	return( pd.DataFrame([geneStarts-1450, geneEnds+1450]).transpose() )


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


	ts = pyslim.load(args.trees)
	N = ts.get_sample_size()

	alive = ts.individuals_alive_at(0)


## Make a dict of the environments
	enviDict = {}
	envi = []
	count = 0
	with open(args.optima) as file:
		for line in file:
			enviDict[ count ] = int( line.strip() )
			count += 1

## Let's make a dict of all the demes and the individuals from each
	subPopDict = {}	
	for a in alive:
		try:
			subPopDict[ ts.individual(a).population ]. append(a)
		except:
			subPopDict[ ts.individual(a).population ] = [a]
		envi.append( enviDict[ts.individual(a).population] )
#	random.shuffle(envi)
#	print(envi)
	
	if len( enviDict.keys() ) != len( subPopDict.keys() ):
		print( "Uneven numbers of subpopulations in the simulation and the optima file that was given" )
		return

	nPops = len( enviDict.keys() )
	
	diploids_per_pop = len( subPopDict[0] )

	print( diploids_per_pop ,"diploid individuals per population")
	print( nPops ,"populations")


## The simulations have a common genome structure, 
## mimicking a group of species with syntenic genomes,
	genome = genomeMaker()
	genome['names'] = np.array(['gene'+str(i) for i in range(genome.shape[0])])

## Make a dict for the environments for the populations in the 
## 10 deme simulations. Make a key for population 99. 
## There are ghost samples from pop 99, so to avoid a KeyError
## I just make a little dict for them.
#	enviDict = {1:-6, 2:-4, 3:-2, 4:-1, 5:0, 6:0, 7:1, 8:2, 9:4, 10:6, 99:99}

## Make a vector that contains the environment for each sample

	variants = {}
		
## Get the allele frequencies for all segregating sites at QTL

	print('extracting segregating sites from trees...')

# Make a little dict that will be populated by the genes that contribute to LA
	adapGenes = {}
	
## Iterate over all variants present on the tree from SLiM

	for vs in ts.variants():
		alleleFreqs = (vs.genotypes[::2] + vs.genotypes[1::2]) ## I hate this line of code...
		# It is doing exactly what I want it to do. 
		if vs.genotypes.sum() == N or vs.genotypes.sum() == 0: continue # Remove fixed mutations

		gene = findGene(vs.position, genome)
		
		if len(vs.site.mutations) >1:
			print('This site has more than a single variant:',vs.position)
		
		md = pyslim.decode_mutation(vs.site.mutations[0].metadata)
		
		selCoeff = md[0].selection_coeff
		
		if sum(alleleFreqs) == N/2:
			continue
		
		adapGenes[gene[0]] = 1 ## This is to keep track of which genes have mutations that affect the phenotype
		
		genos_per_pop = ([i  for i in chunks(alleleFreqs,int(diploids_per_pop/2))])
#		print(genos_per_pop)
		genos = []

		for g in genos_per_pop:
			a1 = (len(g) * 2) - g.sum() # The number of A alleles
			a2 = g.sum() # The number of a alleles
			genos.append( int(a1) )
			genos.append( a2 )
		if sum(genos) != len(alleleFreqs) *2:
#			print(genos_per_pop)
#			print(genos, sum(genos))
			print('something went haywire with the genotype counts')
		variants[vs.position] = [alleleFreqs, selCoeff, genos]
#		print( [alleleFreqs, selCoeff, genos] )

	print(len(variants.keys()), 'segregating sites affecting QTL')		



## Set the mutation rate for neutral mutations 
	mut_rate = 5e-9 # HARD CODED, SO CHANGE IF NECESSARY

## Sprinkle mutations onto the coalescent tree
	print('Sprinkling mutations onto trees')
	sprinkled = msprime.mutate(ts, rate= mut_rate, keep=True)
	

## Iterate over all polymorphisms in the tree and 
	count = 0

	for variant in sprinkled.variants():
		count +=1
		if count % 1000 == 0:
			print('extracted',count,'neutral sites')
		if variant.position in variants.keys():
			s = -99
		else:
			s = 0.0
		alleleFreqs = (variant.genotypes[::2] + variant.genotypes[1::2])
		all_alleles = alleleFreqs.sum()
		if variant.genotypes.sum() == N: continue # Remove fixed mutations
#		print(list(variant.genotypes))
		genos_per_pop = ([i  for i in chunks(alleleFreqs,int(diploids_per_pop/2))])
		
		genos = []

		for g in genos_per_pop:
			a1 = (len(g) * 2) - g.sum()
			a2 = g.sum()
			genos.append( int(a1) )
			genos.append( a2 )

		if sum(genos) != len(alleleFreqs) *2:
			print('something went haywire with the genotype counts')
		
#		print(genos)


#		[(variant.genotypes[i] + variant.genotypes[i+1])/2 for  i in range(N)[::2]]
#		alleleFreqs =  (variant.genotypes[:int(N/2)] + variant.genotypes[int(N/2):] )/2
		variants[variant.position] = [alleleFreqs, s, genos]
	

## For each variant, identify the gene it is within,
## calculate the Spearman's Rho  as well as pq
## These data are then collated and made into a dataFrame
	print(count, "neutral segregating sites" )
	data = []
	count = 0
	print('Performing GEA on the resulting data')
#	output = open(args.output, 'w')
#	output.write('pos,gene,s,rho,rho2,pval,pbar_qbar\n')
	scanners = []
	for v in variants.keys():
		s = variants[v][1]
		alleleFreqs = np.array(variants[v][0])

		gene = findGene(v, genome)

		if not gene:
			continue

		elif gene == -99:
			print(gene)
			print('FAILED')
			return
		
		try:
			LA = adapGenes[gene[0]]
		except KeyError:
			LA = 0	

		pbar = alleleFreqs.sum()/ (len(alleleFreqs)*2)

#		print(alleleFreqs)
		if pbar == 1:continue
		if pbar > 1:
			print("WTF, what is going on here: pbar > 1")
			print(alleleFreqs)
			continue
		pbar_qbar = pbar * ( 1 - pbar )
		maf = min(pbar, 1-pbar)
		if pbar_qbar <0:
			print("WTF, what is going on here: pbar_qbar < 1")
			print(alleleFreqs)
			continue

		spearman = scipy.stats.spearmanr(envi, alleleFreqs)
		k_tau, k_tau_p_value = scipy.stats.kendalltau(envi, alleleFreqs)
		
		if spearman.correlation == 0: continue
			
		dataPoint = {'pos':v,
		'gene':gene[0],
		's':s,
		'spearman_rho':spearman.correlation,
		'spearman_rho2':spearman.correlation**2,
		'spearman_rho_pval':spearman.pvalue,
		'kendall_tau':k_tau,
		'kendall_tau_pval':k_tau_p_value,
		'pbar':pbar,
		'pbar_qbar':pbar_qbar,
		'maf':maf,
		'LA':LA}
		data.append( dataPoint )
		scanners.append(v)
#		if count == 1000:break

## convert the GEA data points into a summary CSV 
	pd.DataFrame(data).sort_values(by=['pos']).to_csv(args.output +'.csv', index = False)
	
## Make a BayScan Input File from the input data

#	BayScan = []
#	for b in sorted(vars.keys()):
#		BayScan.append(vars[b][2])
#		print(b, vars[b][2])
	BayScan = [variants[b][2] for b in sorted(scanners)]
	
	pd.DataFrame(BayScan).to_csv(args.output +'.bayPass.txt', index = False, header = False, sep = ' ')
	
	#.to_csv(args.output, index = False)
main()



