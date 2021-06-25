import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from GEA_functions import *

## Here's a script to calculate the Weighted Z score (and the top-candidate test from the AdapTree data structures - i.e. those that are present in the Dryad repository)

class corLine:

	def __init__(self, c, envIndex):

		dat = c.strip().split()

		self.contig = dat[6]

		self.pos = int(dat[7])

## The data in the Correlations file are stored as contig-pos, so I'll parse that out
## Note that this is VERY sensitive to the file structure staying the same
		if dat[envIndex] == "NA":
				self.pVal = "NA"
		else:
				self.pVal = float(dat[envIndex])

		MAF =  1 - float(dat[13])

		self.pbar_qbar = MAF * (1 - MAF) 

## A generator function to spit out the consecutive lines for the contig
def contigGenerator(correlationFile, envFilter):


	current_contig = ""

	contig_data = []

	with open( correlationFile ) as cor:
		for c in cor:
## Ignore the header
			if c.startswith("X.annotation"):continue
			
			currentLine = corLine(c, envFilter)

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

def contigSnpTable(snps):

	data_for_table = []

	for s in snps:
		
# 		gene_bool = (s.pos>=contig_dat.start)&(s.pos<=contig_dat.end)
# 
# 		if sum(gene_bool) == 0:continue
# 
# 		gene_name =s":".join(list(contig_dat.attribute[gene_bool]))
# 		gene_start = list(contig_dat.start[gene_bool] )[0]
# 
# 		gene_end = list( contig_dat.end[gene_bool] )[0]
# 
# 		if sum(gene_bool) > 1:
# 			print("more than one gene overlapping with SNP:", s.contig, s.pos)

		data_for_table.append({"contig":s.contig,
									"pos":s.pos,
#									"rho":s.rho,
									"pVal":s.pVal,
									"pbar_qbar":s.pbar_qbar,
									"gene":s.contig})

	return pd.DataFrame(data_for_table)	
	
	
def correlationThreshold( corData, targetEnv, percentile_threshold = 99.9):
## Make an empty container to dump the pVals into
	pValues = []
	with open( corData ) as cor:
		for c in cor:
			if c.startswith("X.annotation") or c.startswith("X_annotation"):continue
			currentLine = corLine(c, targetEnv)
			if currentLine.pVal == "NA": continue
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

	parser.add_argument("--output", 

			required = True,

			dest = "output",

			type = str, 

			help = "The name of the output file (the environment will be prepended to the file name so be sure to write to this dir!)")

	parser.add_argument("--env", 

			required = False,

			dest = "env",

			type = str, 

			help = "If you want to analyse just a single environment, give it here [DD_0]",
			
			default = "DD_0")
			
<<<<<<< HEAD
	parser.add_argument("--empirical_p", 

			required = False,

			dest = "empirical_p",

			action = "store_true", 

			help = "[OPTIONAL] Give this flag if you want to analyse empirical p-values that have been appended to the dataframe")
=======
	parser.add_argument("--bay", 

			required = False,

			dest = "bay",

			action = "store_true", 

			help = "[OPTIONAL] Give this flag if the analysis files are BayesFactors from bayEnv.")
>>>>>>> origin/master

	args = parser.parse_args()

## For all environmental variables:
## MAT MWMT MCMT TD MAP MSP AHM SHM DD_0 DDS NFFD bFFP FFP PAS EMT EXT Eref CMD

	adapTreeEnvs = ["LAT"	,"LONG"	,"ELEVATION"	,"MAT"	,"MWMT"	,"MCMT"	,"TD"	,"MAP"	,"MSP"	,"AHM"	,"SHM"	,"DD_0"	,"DD5"	,"NFFD"	,"bFFP"	,"eFFP"	,"FFP"	,"PAS"	,"EMT"	,"EXT"	,"Eref"	,"CMD"	,"Budset_p"	,"Budbreak_p"	,"Height_season_1_p"	,"Height_season_2_p"	,"Diameter_p"	,"Shoot_weight_p"	,"Root_weight_p"	,"Max_growth_rate_p"	,"Linear_growth_days_p"	,"X5_growth_complete_days_p"	,"X95_growth_complete_days_p"	,"X5_95_growth_days_p"	,"Fall_cold_injury_p"	,"Winter_cold_injury_p"	,"Spring_cold_injury_p"	,"root_wt_shoot_wt_p"	,"root_wt_shoot_wt_p_1"]

	adapTreeEnvDict = {}
	
	for i in range(len(adapTreeEnvs)):
		adapTreeEnvDict[adapTreeEnvs[i]] = i + 32 ## to adjust for the length of the adaptree SNP table files

## If you specify just a single environment (which is recommended for parallelisation)
	if args.env:
## Make sure that the environment specified is actually valid
		if args.env not in adapTreeEnvs:
			print("you did not specify a valid environment. Take a look at your command, bozo")
			return
		envs = [args.env]

## We're going to analyse each env. separately
	for env in envs:
<<<<<<< HEAD
## If you've made the file of empirical ps, only analyse the last column
		if args.empirical_p:
=======
		if args.bay:
>>>>>>> origin/master
			envIndex = -1
		else:	
			envIndex = adapTreeEnvDict[env]
		all_contigs = []

## Now let's calculate the outlier threshold (for the TC test) from the data - Make sure it's a percentile!
		threshold_99th = correlationThreshold( args.correlations, envIndex, percentile_threshold = 99)
<<<<<<< HEAD

=======
>>>>>>> origin/master
		if threshold_99th == None:
			print("Something went wrong when identifying the outlier threshold")
			return
		print("99th percentile:",threshold_99th) 


		print("Analysing:",env)
		
		count = 0
		
<<<<<<< HEAD
## Iterate over contigs spat out by the contig generator
		for contig,SNPs in contigGenerator(args.correlations, envIndex):
			count += 1
#			if count ==100: break

## Grab all the genes present on this contig
			contigDF = contigSnpTable(SNPs)
			print(contig, count)
=======
## Iterate over contigs spat out by the 
		for contig,SNPs in contigGenerator(args.correlations, envIndex):
			count += 1
			print(contig, count)
			if count ==100: break

## Grab all the genes present on this contig
			contigDF = contigSnpTable(SNPs)
>>>>>>> origin/master

# Remove all NAs from the SNP set
			contigDF  = contigDF[contigDF["pVal"]!="NA"]
## If there are no annotations on the current contig, move to the next
			if contigDF.shape == (0,0):
				continue
			if contigDF.shape[0] <5: continue
## Get the average position of each annotation - not used for the analysis, just for downstream plotting

			position = contigDF.groupby(["gene"])["pos"].mean().to_frame()

			
## Perform the WZA on the annotations in the contig
			wza = WZA(contigDF, "pVal")
## Perform the top-candidate for the annotations in the contig
			TopCan = top_candidate( contigDF, threshold_99th, 0.01, "pVal", 0.9999, MAF_filter = 0.05 )
## Combine the results
			result = pd.concat([ position, wza, TopCan] , axis = 1, sort = True ).reset_index()

			result["contig"] = contig

			result["env"] = env

			all_contigs.append( result )

			if count%1000 == 0:
				print( count,"contigs analysed")
			
			
		if len( all_contigs ) == 0: continue 
		
## Combine all contig-specific dataframes into a single big one
		outputDF = pd.concat(all_contigs)

## Get the expected proportion of hits per hit-bearing gene
		expected = outputDF[outputDF["hits"]!=0]["hits"].sum() / outputDF[outputDF["hits"]!=0]["SNPs"].sum()

		print("expected proportion", expected)
		top_candidate_p = [ scipy.stats.binom_test(h, s, expected, alternative = "greater" ) for h, s in zip(outputDF.hits, outputDF.SNPs)]
		
		outputDF["top_candidate_p"] = top_candidate_p 
#		outputDF["top_candidate_p"] = scipy.stats.binom_test(outputDF.hits, outputDF.SNPs, expected, alternative = "greater" ) 
		expected_hits = [ scipy.stats.binom.ppf( 0.9999 , s, expected )  for s in outputDF.SNPs]

		outputDF["expected_hits"] =  scipy.stats.binom.ppf( 0.9999 , outputDF.SNPs, expected )  

## Calculate the top-candidate index
#		outputDF = pd.concat(all_contigs)		


## Write the dataframe to an output file
		outputDF.to_csv(env + "_" + args.output,index = False)


main()




