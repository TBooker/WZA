import msprime, pyslim, random, allel
import pandas as pd
import scipy.stats
from scipy.spatial import distance_matrix
import numpy as np
import argparse, tskit
import matplotlib.pyplot as plt



def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--trees", 
			required = True,
			dest = "trees",
			type = str, 
			help = "Coalescent trees for the simulation")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "What name do you want to give to the output file? [DON'T GIVE FILE EXTENSION, THE SCRIPT MAKES TWO OUTPUT FILES")
	parser.add_argument("--island", 
			dest = "island",
			action = "store_true", 
			help = "Use this flag if the data comes from the Island model")
	args = parser.parse_args()

	wins = [i for i in range(0, 1000000+3000, 1000)] 
	ts = pyslim.load(args.trees)

#	print("recap")

	mut_rate = 1e-7 # HARD CODED, SO CHANGE IF NECESSARY

	pops = {} 
	for i in ts.individuals_alive_at(0):
		if ts.individual(i).population == 999: continue

		try: 
			pops[ts.individual(i).population].append( i )
		except KeyError:
			pops[ts.individual(i).population] = [i]
	
# Let's take a sample of 40 pops and 20 diploids from each and make a new tree from them
	sampled_indivs = []
	for i in np.random.choice( list( pops.keys() ), 1, replace = False):
		for j in np.random.choice( np.array(pops[i]), 10 , replace = False):
			sampled_indivs.extend( ts.individual(j).nodes )
			
	r2_tree2 = ts.simplify(sampled_indivs)

		#print(i.population)
#	print(ts.sample_size)
## Sprinkle mutations onto the coalescent tree
#	print('Sprinkling mutations onto trees')

#	print("sprinkle 1")
	sprinkled = msprime.mutate(r2_tree2, rate= mut_rate, keep=True)
#	MAF_005 = []  
#	for i in sprinkled.variants():
#		print(i.position)
#		if sum(i.genotypes)/len(i.genotypes) > 0.01:
#			MAF_005.append(1)
#		else:
#			MAF_005.append(0)
#	print( np.array(MAF_005).sum() )		
	

#	return
#	print("sprinkle 2")

#	sprinkled_recap = msprime.mutate(ts, rate= mut_rate, keep=True)

####	print(sprinkled.diversity( windows = 	[0,1e6,1e6+3000]) )

#	sprinkled_recap.diversity( windows = 	wins+ [1003000]
#))
#	plt.step( np.array(wins), sprinkled.diversity( windows = 	wins+ [1003000top
#]) , "xkcd:crimson" ) 
#	plt.step( np.array(wins), sprinkled_recap.diversity( windows = 	wins+ [1003000]) , "peachpuff" ) 
#	plt.show()

	
#	return
	positions = [i.position for i in sprinkled.variants() if i.position < 1e4]
#	print(len(positions))

	LD_interval = sprinkled.keep_intervals(np.array([[0,20000]])) 
	LD = tskit.LdCalculator(LD_interval)
	LD_mat = LD.r2_matrix()
	


	for i in range(len(positions)):
		for j in range(len(positions)):
			if j>=i:
				continue
			else:
				pw_distances_ij = abs( positions[i] - positions[j] )
				if pw_distances_ij > 100000: continue
				print(pw_distances_ij, LD_mat[i,j])
	
main()
