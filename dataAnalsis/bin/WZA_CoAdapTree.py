import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from GEA_functions import *

## Here's a script to calculate the Weighted Z score (and the top-candidate test from the CoAdapTree data structures)

class corLine:

	def __init__(self, c):

		dat = c.strip().split()

		self.contig = dat[0]

		self.pos = int(dat[1])

## The data in the Correlations file are stored as contig-pos, so I'll parse that out
## Note that this is VERY sensitive to the file structure staying the same
		self.env = dat[2]

		self.rho = float(dat[3])

		self.pVal = float(dat[4])

		MAF = float(dat[5])

		self.pbar_qbar = MAF * (1 - MAF) 

		self.emp_pVal = float(dat[6])


## A generator function to spit out the consecutive lines for the contig
def contigGenerator(correlationFile, envFilter):

	current_contig = ""

	contig_data = []

	with open( correlationFile ) as cor:
		for c in cor:
## Ignore the header
			if c.startswith("snp") or c.startswith("contig"):continue
			
			currentLine = corLine(c)

			if currentLine.env != envFilter: continue

			if currentLine.contig != current_contig:

				if current_contig == "":

					current_contig = currentLine.contig

					contig_data.append( currentLine )

				else:

					yield current_contig, contig_data

					contig_data = [currentLine]

					current_contig = currentLine.contig

			elif currentLine.contig == current_contig:

				contig_data.append( currentLine )

	yield current_contig, contig_data

def contigSnpTable(snps, contig_dat):

	data_for_table = []

	for s in snps:

		gene_bool = (s.pos>=contig_dat.start)&(s.pos<=contig_dat.end)

		if sum(gene_bool) == 0:continue

		gene_name = ":".join(list(contig_dat.attribute[gene_bool]))

		gene_start = list(contig_dat.start[gene_bool] )[0]

		gene_end = list( contig_dat.end[gene_bool] )[0]

		if sum(gene_bool) > 1:
			print("more than one gene overlapping with SNP:", s.contig, s.pos)

		data_for_table.append({"contig":s.contig,
									"pos":s.pos,
									"rho":s.rho,
									"pVal":s.pVal,
									"empirical_pVal":s.emp_pVal,
									"pbar_qbar":s.pbar_qbar,
									"gene":gene_name,
									"gene_start":gene_start,
									"gene_end":gene_end})

	return pd.DataFrame(data_for_table)	
	
	
def correlationThreshold( corData, targetEnv, percentile_threshold = 99.9):
## Make an empty container to dump the pVals into
	pValues = []
	with open( corData ) as cor:
		for c in cor:
			if c.startswith("snp"):continue
			currentLine = corLine(c)
			if currentLine.env != targetEnv: continue
			else:
				pValues.append(currentLine.pVal)
	if len(pValues) == 0:
		return None
	else:
		return( 1 -np.percentile( 1 - np.array(pValues), percentile_threshold ) )


def main():

## Define command line args

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--correlations", "-c", 

			required = True,

			dest = "correlations",

			type = str, 

			help = "The file containing the correlations")

	parser.add_argument("--annotations","-a", 

			required = True,

			dest = "annotations", 

			help = "The file of annotations. The script asssumes GFF as default, set the '--bed' flag if using that format")

	parser.add_argument("--output", 

			required = True,

			dest = "output",

			type = str, 

			help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

	parser.add_argument("--env", 

			required = False,

			dest = "env",

			type = str, 

			help = "[OPTIONAL] If you want to analyse just a single environment, give it here")

	parser.add_argument("--bed", 

			required = False,

			dest = "bed",

			action = "store_true", 

			help = "[OPTIONAL] Give this flag if the analysis files are in BED format. Otherwise the script assumes GFF format")

	args = parser.parse_args()
	
	if args.bed:
		annotations = pd.read_csv(args.annotations , 

					sep = "\t",

						names = ["seqname",

								"start",

								"end",

								"attribute"])
## Add 1 to the positions to make correct for 0-based BedTools
		annotations["start"] +=1

	else:
## GFF header from ENSEMBL webpage
		annotations = pd.read_csv(args.annotations , 

					sep = "\t",

						names = ["seqname",

								"source",

								"feature",

								"start",

								"end",

								"score",

								"strand",

								"frame",

								"attribute"])

## For all environmental variables:
## MAT MWMT MCMT TD MAP MSP AHM SHM DD_0 DDS NFFD bFFP FFP PAS EMT EXT Eref CMD

	envs =["MAT","MWMT","MCMT","TD",
			"MAP","MSP","AHM","SHM",
			"DD_0","DD5","NFFD","bFFP",
			"FFP","PAS","EMT","EXT","Eref","CMD"]
## If you specify just a single environment (which is recommended for parallelisation)
	if args.env:
## Make sure that the environment specified is actually valid
		if args.env not in envs:
			print("you did not specify a valid environment. Take a look at your command, bozo")
			return
		envs = [args.env]

## We're going to analyse each env. separately
	for env in envs:
		all_contigs = []

## Now let's calculate the outlier threshold (for the TC test) from the data - Make sure it's a percentile!
		threshold_99th = correlationThreshold( args.correlations, env, percentile_threshold = 99)
		if threshold_99th == None:
			print("Something went wrong when identifying the outlier threshold")
			return
		print("99th percentile:",threshold_99th) 


		print("Analysing:",env)
		
## Iterate over contigs spat out by the 
		for contig,SNPs in contigGenerator(args.correlations, env):
## Grab all the genes present on this contig
			contigDF = contigSnpTable(SNPs, annotations[annotations["seqname"] == contig])
		#	print(contigDF)
## If there are no annotations on the current contig, move to the next
			if contigDF.shape == (0,0):
				continue

## Get the average position of each annotation - not used for the analysis, just for downstream plotting

			position = contigDF.groupby(["gene"])["pos"].mean().to_frame()

			anno_start = contigDF.groupby(["gene"])["gene_start"].mean().to_frame()
			
			anno_end = contigDF.groupby(["gene"])["gene_end"].mean().to_frame()


## Perform the WZA on the annotations in the contig using parametric ps
			wza_pVal = WZA(contigDF, "pVal", varName = "Z_p")

## Perform the WZA on the annotations in the contig using empirical ps
			wza_emp_pVal = WZA(contigDF, "empirical_pVal", varName = "Z_empP")

## Perform the top-candidate for the annotations in the contig
			TopCan = top_candidate( contigDF, threshold_99th, 0.01, "pVal", 0.9999, MAF_filter = 0.05 )

## Combine the results
			result = pd.concat([ position, wza_pVal, wza_emp_pVal, TopCan, anno_start, anno_end] , axis = 1, sort = True ).reset_index()

			result["contig"] = contig

			result["env"] = env

			all_contigs.append( result )
			
		if len( all_contigs ) == 0: continue 
		
## Combine all contig-specific dataframes into a single big one
		outputDF = pd.concat(all_contigs)
		
## Get the expected proportion of hits per hit-bearing gene
		expected = outputDF[outputDF["hits"]!=0]["hits"].sum() / outputDF[outputDF["hits"]!=0]["SNPs"].sum()

		print("expected proportion", expected)

		top_candidate_p = [ scipy.stats.binom_test(h, s, expected, alternative = "greater" ) for h, s in zip(outputDF.hits, outputDF.SNPs)]
		
		outputDF["top_candidate_p"] = top_candidate_p 

#		expected_hits = [ scipy.stats.binom.ppf( 0.9999 , s, expected )  for s in outputDF.SNPs]

		outputDF["expected_hits"] =  scipy.stats.binom.ppf( 0.9999 , outputDF.SNPs, expected )  

		
## Write the dataframe to an output file
		outputDF.to_csv(env + "_" + args.output,index = False)


main()
