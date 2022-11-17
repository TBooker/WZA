import msprime, pyslim, random, allel, glob
import pandas as pd
import scipy.stats
from scipy.spatial import distance_matrix
import numpy as np
import argparse, tskit
import matplotlib.pyplot as plt


def getLD(tree_sequence, output_name, selection = False):
	if selection == True:
		sel_pos_raw = []
		for v in tree_sequence.variants():
			if v.genotypes.sum() < 4000: continue ## A rough way of removing low frueqnecy alleles
			else: 
				sel_pos_raw.append(round(v.position/10000)*10000)
		LD_intervals = [ [s, s+9999] for s in list(set(sel_pos_raw))]

	else:
		## these genes always evolve neutrally
		LD_intervals = [[9000000,9010000-1], [9200000,9210000-1], [9400000,9410000-1], [9600000,9610000-1], [9800000,9810000-1]]

	print(LD_intervals)
	pops = {} 
	for i in tree_sequence.individuals_alive_at(0):
		if tree_sequence.individual(i).population == 999: continue

		try: 
			pops[tree_sequence.individual(i).population].append( i )
		except KeyError:
			pops[tree_sequence.individual(i).population] = [i]
	
# Let's take a sample of 10 individuals from 1 pops and make a new tree from them
	sampled_indivs = []
	for i in np.random.choice( list( pops.keys() ), 1, replace = False):
		for j in np.random.choice( np.array(pops[i]), 10 , replace = False):
			sampled_indivs.extend( tree_sequence.individual(j).nodes )
	
	
	r2_tree2 = tree_sequence.simplify(sampled_indivs)
	mut_tree = msprime.mutate(r2_tree2, 1e-7)

#	print(pops)
#	print(ts.sample_size)
## Sprinkle mutations onto the coalescent tree
#	print('Sprinkling mutations onto trees')

#	print("sprinkle 1")
#	sprinkled = msprime.mutate(r2_tree2, rate= mut_rate, keep=True)



	LD_output = open(output_name, "w")
	for inter in range(len(LD_intervals)):
		interval = LD_intervals[inter]
		positions = [i.position for i in mut_tree.variants() if i.position >= interval[0] and i.position <= interval[1]]
		LD_interval = mut_tree.keep_intervals(np.array([ interval ])) 
		LD = tskit.LdCalculator(LD_interval)
		LD_mat = LD.r2_matrix()
		


		for i in range(len(positions)):
			for j in range(len(positions)):
				if j>=i:
					continue
				else:
					pw_distances_ij = abs( positions[i] - positions[j] )
					if pw_distances_ij > 10000: continue
					LD_output.write(",".join([str(inter), str(pw_distances_ij), str(LD_mat[i,j]) ] ) + "\n")
	LD_output.close()
	

def summary_stats(tree_sequence_file, outputFileName, selection = False):
	ts = pyslim.load(tree_sequence_file)
	getLD(ts, outputFileName, selection )
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
	parser.add_argument("--selection", 
			required = False,
			dest = "selection",
			action = "store_true",
			help = "Do you want to analyse selected sites?")
	args = parser.parse_args()
	
	count = 0
	for i in glob.glob(args.trees + "/*trees"):
		if count >= args.num:
			print("we have done it!")
			break
		name = i.split("/")[-1].split(".")
		direc = "/".join(i.split("/")[0:-1])
		rep = name[0]		
		gen = name[1]
		if args.gen:
			if int(gen) != args.gen:
				continue
		else:
			pass
		print(i)
			
		if args.selection:
			outputFileName = "selection_" + rep + "." + gen + ".LD.csv"

		else:
			outputFileName = rep + "." + gen + ".LD.csv"

		summStats = summary_stats(i, outputFileName, selection = args.selection)
		count += 1

main()
