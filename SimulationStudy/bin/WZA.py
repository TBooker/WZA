import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse

## Pooja, just so it's clear. I normally try and do 2 comments (e.g. ##) for code annotation and reserve 1 (e.g. #) for commenting out code

## I'm going to write a set of functions that perform the statistical tests on each of the genes - looking for signatures of local adaptation
## Each one is going to return a value for each gene - whether a particular test identified positive selection or not. 

## Here's the function that implements the Weighted Z Analysis (WZA) on the SLiMulated data


def WZA( gea , statistic , MAF_filter = 0.05, varName = "Z"):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)
	gea[statistic] = gea[statistic].clip(lower = 1e-15 )
	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

## convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z"] = scipy.stats.norm.ppf(1 - gea[statistic])

## Apply the MAF filter
	gea_filt = gea[ gea["maf"] > MAF_filter ].copy()

## Calculate the numerator and the denominator for the WZA
	gea_filt[varName+"weiZ_num"] = gea_filt["pbar_q_bar"] * gea_filt["z"]
	gea_filt[varName+"weiZ_den"] =  gea_filt["pbar_q_bar"]**2 	
	numerator = gea_filt.groupby(["gene"])[varName+"weiZ_num"].sum().to_frame()
	denominator = np.sqrt(gea_filt.groupby(["gene"])[varName+"weiZ_den"].sum()).to_frame()

## We've calculated the num. and the den., let's make a dataframe that has both 
	weiZ  = pd.concat([numerator,denominator], axis = 1, sort = False)
## Actually calculate the Z scores for each gene
	weiZ[varName] = weiZ[varName+"weiZ_num"] / weiZ[varName+"weiZ_den"]

## One might be interested in calculating a p_value from the Z-scores (though this only works if the data are normal, which they won't be if there's populations structure).
	weiZ[varName+"_hits"] = (weiZ[varName] > scipy.stats.norm.ppf(1 - 0.05/50)).astype(int)

## Return the final dataframe
	return weiZ

## Generic function to calculate a threshold from a set of GEA files and the p-values within
def getThresh( gea_files , percentile, statistic= "pop_k_tau_p_value" ):
## Create an empty list to dump all the p_values into 
	p_vals = []
	for i in gea_files:
# Read in the GEA data 
		gea_data = pd.read_csv( i )

## Convert the p_values into -log10(p_values)
#		gea_data["temp"] = -np.log10( gea_data[statistic])
		p_vals += list( -np.log10( gea_data[statistic] ) )

#		for i in list( gea_data["pop_k_tau_p_value"]  ):
#			print(i , -np.log10(i) )
#			if i == -0:
#				print("!!!!!!!!!!!!!!!")
#				break

#	print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#	print( np.nanpercentile( p_vals, percentile ) )
#	print( np.percentile( p_vals, percentile ) )

	return np.nanpercentile( p_vals, percentile ) 

## Generic function calculate the proportion of values with a a p_value passing a given threshold

def hitsFraction( gea_files, thresh):
	nHits = 0
	nSNPs = 0
	for i in gea_files:
		gea_data = pd.read_csv( i )
		hits =  -np.log10(gea_data["pop_k_tau_p_value"]) > thresh 
		nHits += sum(hits)
		nSNPs += gea_data.shape[0] 
	## Calculating the overall number of hits, or the number of hits per simluation replicate give vrey similar answers. 
	## I'll just go with the overall number for the time being
	return nHits / nSNPs

## A really straightforward GEA classification: Is there a SNP with a small p_value in a given gene?
def simple_classifier( gea, thresh, MAF_filter = 0.05):

## First convert p_values to -log10(p_values), then test against threshold then convert to 1s and 0s
	gea["hits"] = ( -np.log10(gea["pop_k_tau_p_value"]) > thresh).astype(int)
## Drop any SNP that does not meet the MAF filter
	gea_filt = gea[ gea["maf"] > MAF_filter ]

## Count the hits per gene
	num_hits = gea_filt.groupby(["gene"])["hits"].sum().to_frame()
## Make a dataframe for the big hits
	BH  = pd.concat([num_hits], axis = 1, sort = False) 
	BH.columns = ["bigHits"]
	return(BH)

## A function for performing the top-candidate test

def top_candidate( gea, thresh, prop_hits, MAF_filter = 0.05):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit

## Identifty the hits
	gea["hits"] = ( -np.log10(gea["pop_k_tau_p_value"]) > thresh).astype(int)

## Apply the MAF filter
	gea_filt = gea[ gea["maf"] > MAF_filter ]
## Count the hits per gene
	num_hits = gea_filt.groupby(["gene"])["hits"].sum().to_frame()
## Count the SNPs per gene
	num_SNPs = gea_filt.groupby(["gene"])["hits"].count().to_frame()

## Make a single DF with the hits and the SNPs
	TC  = pd.concat([num_hits, num_SNPs], axis = 1, sort = False) 
## Name the cols
	TC.columns = ["hits", "SNPs"]
## Init an empty vector for p_values (the top-candidate index)
	p_vals = []
	for index, row in TC.iterrows():
		p_vals.append( scipy.stats.binom_test(row.hits, row.SNPs, prop_hits, alternative = "greater" ) )
## Add the TC p_vals to the genes
	TC["top_candidate_p"] = p_vals
## Return the resulting dataFrame
	return(TC)

## For each gene, figure out the minimum distance to the nearest gene that explains a substntial portion of phenotypic variance
def distance_to_PVE(localAdaptation, geneMiddles, PVE_threshold):
## Get the gene positions of those that contribute substnatially to phenotypic variation
	PVE_genes_positions = []
	for k in localAdaptation.keys():
		if localAdaptation[k] >= PVE_threshold:
			PVE_genes_positions.append( geneMiddles[k] )
	
## Init a dict to store the gene positions
#	print(PVE_genes_positions)
	dist_dict = {}
## If there is no PVE genes in the simulation, say that the dist is 1e6 (i.e. the length of the chromosome)
	if len(PVE_genes_positions) == 0:
		for k in geneMiddles.keys():
			dist_dict[k] = 1e6 
## If there IS PVE genes in the simulation, get the distance to the nearest one for each gene
	else:
		for k in geneMiddles.keys():
			dist_dict[k] = min([abs(geneMiddles[k] - X) for X in PVE_genes_positions] ) 
## Turn the distance dict into a dataframe
	return( pd.Series( dist_dict , name = "distToPVE") )

## A master function to analyse a directory full of simulation results
def analyseSimSet(simulationSet, fileSet, PVE_threshold = 0.005, bayPass = False):

## Get the 1 in 100 threshold
	threshold_99 = getThresh( simulationSet ,99)
## Get the 1 in 1000 threshold
	threshold_999 = getThresh( simulationSet ,99.9)

## Determine the average fraction of hits per gene
	fracking_hits = hitsFraction( simulationSet  , threshold_99)

## Iterate over all GEA files in the directory
	for r in simulationSet:
		print(r)
		gea_data = pd.read_csv( r )
## IF there are BayPass results in the directory, this is how we handle them...
		if bayPass:
			bayPass_file = r.split(".csv")[0] + "_summary_betai_reg.out"

## BayPass output files are written with a variable number of white spaces as a delimiter
## The regex r"\s*" tells pandas that there are a variable number of white spaces

## The TRY/EXCEPT conditional here tests if a BayPass file is present corresponding to a particular simulation. If not, don't exit, just go to the next simulation
			try:
				bayPass_data = pd.read_csv( bayPass_file ,sep = r"\s*")
			except FileNotFoundError:
				continue

## If we have BayPass results, make a strimmed down version of the data file
## BayPass output header:
## ["COVARIABLE", "MRK", "M_Pearson" ,"SD_Pearson", "BF(dB)", "Beta_is", "SD_Beta_is", "eBPis"]
## If the number of rows inthe GEA file and the BayPass file don't match then that's a problem 
			bayPass_data["pVal"] = 10**-bayPass_data["eBPis"]
			bayPass_data["gene"] = gea_data["gene"]
			bayPass_data["maf"] = gea_data["maf"]
			bayPass_data["pbar_q_bar"] = gea_data["pbar_q_bar"]

## Make a dataframe with all the PVEs per gene
		LA = gea_data.groupby(["gene"])["LA"].mean().to_frame()
## Make a dataframe with all the positions per gene (average position across SNPs)
		position = gea_data.groupby(["gene"])["position"].mean().to_frame()
## Make a dataframe with all the distances to the nearest PVE gene
		distance = distance_to_PVE(LA.to_dict()["LA"], position.to_dict()["position"], PVE_threshold)

## Do the top-candidate and WZA on uncorrected data
## Do WZA on BayPass - I've not yet implemented top-candidate for the BayPass results. I'll do this soon, but hopefully that's not a stumbling block for the time being.
		if bayPass:
			wza_results_kendall = WZA( gea_data ,  "pop_k_tau_p_value", MAF_filter = 0.05, varName = "Z_kendall") 
			wza_results_bayPass = WZA( bayPass_data , "pVal", MAF_filter = 0.05, varName = "Z_bayPass") 
		else:
			wza_results = WZA( gea_data ,  "pop_k_tau_p_value", MAF_filter = 0.05) 

		TC_results = top_candidate( gea_data, threshold_99, fracking_hits ) 

		simple_results = simple_classifier( gea_data, threshold_999 )

		if bayPass:
			results = pd.concat([ distance, wza_results_kendall, wza_results_bayPass, TC_results, LA, position, simple_results ] , axis = 1, sort = True ).reset_index()
		else:
			results = pd.concat([ distance, wza_results, TC_results, LA, position, simple_results ] , axis = 1, sort = True ).reset_index()
		results["rep"] = fileSet[r]["rep"]
		results["sig_a"] = fileSet[r]["sig_a"]
		results["U_a"] = fileSet[r]["U_a"]
		results["Vs"] = fileSet[r]["Vs"]
#		all_results.append( results )
		yield( results )

""" R Code for the WZA would look something like:
gea$z <- qnorm(gea$pop_k_tau_p_value, lower.tail = F)
gea$z[gea$z==-Inf] <- qnorm(0.999, lower.tail = F)

gea_filt<-gea#[gea$maf > 0.01, ]

gene_pos <- tapply( gea_filt$position, as.factor(gea_filt$gene), mean)

LA <- tapply( gea_filt$LA, as.factor(gea_filt$gene), mean)

weiZ_num <- tapply( gea_filt$pbar_q_bar * gea_filt$z, as.factor(gea_filt$gene), sum)

weiZ_den <- sqrt( tapply( gea_filt$pbar_q_bar**2, gea_filt$gene, sum))

Z <- weiZ_num/weiZ_den
"""

def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--csv", 
			required = True,
			dest = "csv",
			type = str, 
			help = "The directory containing the CSVs")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "The name of the output file")
	parser.add_argument("--neutral", 
			action = "store_true",
			dest = "neutral", 
			help = "Were these neutral simulations?")
	parser.add_argument("--sampled", 
			action = "store_true",
			dest = "sampled", 
			help = "Were these simulations downsampled (i.e. is there the tag 'd40_n50' on the file names?)")
	parser.add_argument("--bayPass", 
			action = "store_true",
			dest = "bayPass", 
			help = "Do you want to niclude BayPass results?")			
	args = parser.parse_args()

## Make a dict for the GEA files
	files = {}
## Make a dict for each of the simulation parameter sets
	simSets = {}
## Were the GEA data sampled or do they respresent all demes?
	if args.sampled:
		searchString = "/*_d40_n50.csv"
	else:
		searchString = "/*00.csv"
## What follows here is a lot of string parsing games that I play to extract metadata from file names		
	for i in glob.glob(args.csv + searchString):
		if args.neutral:
			nameString = i.split(".csv")[0].split("/")[-1] 
			rep = nameString.split("_")[0]
			sig_a = "0"
			U_a = "0"
			Vs = "0"
		else:
			nameString = i.split(".csv")[0].split("/")[-1] 
			rep = nameString.split("_")[0].split(".")[0]
#			sig_a = ".".join( nameString.split("_")[0].split(".")[1:] )
			sig_a = nameString.split("_")[1]
			U_a = nameString.split("_")[2]
#			Vs = args.Vs
			Vs = nameString.split("_")[3]
#			print(nameString, rep,sig_a, U_a, Vs)

		files[i] = {"rep":rep, "sig_a":sig_a, "U_a":U_a, "Vs":Vs }

		try:
			simSets[sig_a + "_" + U_a + "_" + Vs].append(i)
		except KeyError:
			simSets[sig_a + "_" + U_a + "_" + Vs] = [i]

	all_results = []

	for s in simSets.keys():
		for q in analyseSimSet(simSets[s], files, bayPass = args.bayPass):
			all_results.append( q )
			#break
		#break
#		all_results += analyseSimSet(simulationSet)


	myData = pd.concat( all_results, sort = True ).reset_index(drop=True)
	myData.to_csv( args.output, index = False)


main()
