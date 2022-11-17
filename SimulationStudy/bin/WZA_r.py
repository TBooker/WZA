import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from collections import Counter

## I'm going to write a set of functions that perform the statistical tests on each of the genes - looking for signatures of local adaptation
## Each one is going to return a value for each gene - whether a particular test identified positive selection or not. 

## Here's the function that implements the Weighted Z Analysis (WZA) on the SLiMulated data


def WZA( GEA , statistic , MAF_filter = 0.05, varName = "", groupVar = "gene"):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	GEA[statistic] = GEA[statistic].clip(lower = 1e-15 )
	GEA[statistic] = GEA[statistic].replace(1 , 1-1e-3)

## convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	GEA["z"] = scipy.stats.norm.ppf(1 - GEA[statistic])

## Apply the MAF filter
	gea_filt = GEA[ GEA["maf"] > MAF_filter ].copy()

## Calculate the numerator and the denominator for the WZA
	gea_filt[varName+"_weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z"]
	gea_filt[varName+"_weiZ_den"] =  gea_filt["pbar_qbar"]**2 	
	numerator = gea_filt.groupby([groupVar])[varName+"_weiZ_num"].sum().to_frame()
	denominator = np.sqrt(gea_filt.groupby([groupVar])[varName+"_weiZ_den"].sum()).to_frame()

## We've calculated the num. and the den., let's make a dataframe that has both 
	weiZ  = pd.concat([numerator,denominator], axis = 1, sort = False)
## Actually calculate the Z scores for each gene
	weiZ[varName + "_Z"] = weiZ[varName+"_weiZ_num"] / weiZ[varName+"_weiZ_den"]

## Return the final dataframe
	return weiZ
## Generic function calculate the proportion of values with a a p_value passing a given threshold


## A really straightforward per SNP classification: 
## What is the rank order of genes based on individual SNP values?

## R Code
"""

NumberOfTruePositives <- function( gea_results, path, suffix, LA_thresh ){
  falsePositiveDF <- data.frame(0,0,0,0,0)
  names(falsePositiveDF) <-  c("rep","top", "truePositive_WZA","truePositive_TC","truePositive_SNPs")
  
  for (rep in 1:20){
    
    temp <- gea_results[gea_results$rep == rep,]
    temp <- temp[temp$LA != -99,]
    numberLA <- length(droplevels(temp[temp$LA > LA_thresh,]$gene))
    LA_genes <- droplevels(unique(temp[temp$LA > LA_thresh,]$gene))
    
    snps <- read.csv(paste(path, rep, suffix, sep = ''))
    snps <- snps[snps$LA != -99,]
    
    snps <- snps[ with(snps, order(ave(geno_k_tau_p_value, gene, FUN = min), geno_k_tau_p_value)),]
    snp_based_genes <- unique( snps$gene )
    
    for (n in 1:50){
      
      WZA_slice <- temp[temp$Z_empirical_p <= (n+1)/1000,]
      TC_slice <- temp[temp$TC_empirical_p <= (n+1)/1000,]
      SNP_slice <- snp_based_genes[1:n]
      
      #      SNP_slice <- snps[snps$rank <= n,] ## Grab the top n SNPs from the genome scan
      
      #      numberLA_SNP_slice <- length(unname(unique(SNP_slice[-log10(SNP_slice$LA) > 3.0,]$gene)))
      
      truePositive_WZA = sum( WZA_slice$LA > LA_thresh )/numberLA
      truePositive_TC = sum( TC_slice$LA > LA_thresh )/numberLA
      truePositive_SNPs = sum(SNP_slice%in%LA_genes)/numberLA
      
      falsePositiveDF[(50*(rep-1))+n,] <- c(rep, n, truePositive_WZA, truePositive_TC, truePositive_SNPs)
      
    }
  }
  falsePositiveDF
}

"""

def simple_classifier( gea, statistic, MAF_filter = 0.05, bayPass = False, label = "", groupVar = "gene"):

	gea_filt = gea[ gea["maf"] > MAF_filter ]
	
	## Get the maximum test statistic per gene
	## A little Pandas magic...
		# Group by "gene"
		# Get the maximum test statistic per gene
	genes_ranked = pd.DataFrame(  gea_filt.groupby([groupVar])[statistic].min() )

	genes_ranked = genes_ranked.rename( columns = {statistic : label+"_"+statistic} )

	genes_ranked[label+"_rank"] = genes_ranked[label+"_"+statistic].rank()

	return genes_ranked	

## A function for performing the top-candidate test
def top_candidate( gea_raw, statistic, MAF_filter = 0.05, log10= True, label = ""):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probability of getting a hit
	

## Apply the MAF filter
	gea = gea_raw[gea_raw["maf"] > MAF_filter ].copy()


	if log10:
		thresh = np.quantile(-np.log10(gea[statistic]),0.99)
		gea["hits"] = ( -np.log10(gea[statistic]) > thresh).astype(int)
#		gea["hits"] = ( -np.log10(gea[statistic]) > np.quantile(-np.log10(gea[statistic]),0.99)).astype(int)
	else:
		thresh = np.quantile(gea[statistic],0.99)
		gea["hits"] = ( gea[statistic] > thresh).astype(int)
	prop_hits = 0.01

	
## Count the hits per gene
	num_hits = gea.groupby(["gene"])["hits"].sum().to_frame()
## Count the SNPs per gene
	num_SNPs = gea.groupby(["gene"])["hits"].count().to_frame()
## Make a single DF with the hits and the SNPs
	TC  = pd.concat([num_hits, num_SNPs], axis = 1, sort = False) 
## Name the cols
	if label != "":
		name_to_add = label +"_"
	else:
		name_to_add = ""

	TC.columns = [name_to_add + "hits",  name_to_add +"SNPs"]

	aveHits = TC[TC[name_to_add +"hits"]!=0][name_to_add +"hits"].sum() / TC[TC[name_to_add +"hits"]!=0][name_to_add +"SNPs"].sum() 

## Init an empty vector for p_values (the top-candidate index)
	p_vals = []
	p_vals_TC2 = []

	for index, row in TC.iterrows():
		p_vals.append( scipy.stats.binom_test(row[name_to_add + "hits"], row[name_to_add + "SNPs"], prop_hits, alternative = "greater" ) )
		p_vals_TC2.append( scipy.stats.binom_test(row[name_to_add + "hits"], row[name_to_add + "SNPs"], aveHits, alternative = "greater" ) )

## Add the TC p_vals to the genes
	TC[name_to_add + "top_candidate_p"] = p_vals
	TC[name_to_add + "top_candidate_p_TC2"] = p_vals_TC2

## Return the resulting dataFrame
	return(TC)


## For each gene, figure out the minimum distance to the nearest gene that explains a substntial portion of phenotypic variance
def distance_to_PVE(localAdaptation, geneMiddles, PVE_threshold, directional = True):
## Get the gene positions of those that contribute substnatially to phenotypic variation
#	print(localAdaptation, PVE_threshold)
	PVE_genes_positions = []
	for k in localAdaptation.keys():
		if directional:
			if localAdaptation[k] != 0 and localAdaptation[k] >= 1:
				PVE_genes_positions.append( geneMiddles[k] )
		else:
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
def analyseSimSet(simulationSet, PVE_threshold = 0.001, bayPass = "none", IBA = False):

## Iterate over all GEA files in the directory
	for r in simulationSet:
		rep = r.split("/")[-1].split("_")[0]
		print("Analysing replicate:",rep)
## Read in the uncorrected GEA

		if IBA:
			gea_data_raw = pd.read_csv( r )
## If analysing Isolation-by-Adaptation:
## Need to extract those sites that correspond to the neutral chromosome and reset the index of the dataframe
			gea_data = gea_data_raw[gea_data_raw["position"] > 8000000].copy().reset_index()
		else:
			gea_data = pd.read_csv( r )

## IF there are BayPass results in the directory, this is where we handle them...
		if bayPass != "none":
			bayPass_file = bayPass + "/" + r.split("/")[-1].split(".csv")[0] + ".bayPass_summary_betai_reg.out"
			if IBA:
				bayPass_file = bayPass + "/" + rep + "_0.003_directiona_d40n50_IBA_i10000.neutralChrom.bayPass_summary_betai_reg.out"

## BayPass output files are written with a variable number of white spaces as a delimiter
## The regex r"\s+" tells pandas that there are a variable number of white spaces
			bayPass_data = pd.read_csv( bayPass_file ,sep = r"\s+", engine = "python")

## If the number of rows inthe GEA file and the BayPass file don't match then that's a problem 
			if  gea_data.shape[0] !=  bayPass_data.shape[0] :
				print("Something went haywire with the BayPass and Uncorrected data")
				return
				

## Convert the Bayes Factors into empirical P values
## I just rank the values using the rankdata() function (which does so in ascending order) and divide the ranks by the number of rows in the dataframe
## one minus that value gives the empirical P value
			bayPass_data["pVal"] =  1 - scipy.stats.rankdata( bayPass_data["BF(dB)"] )/ len(bayPass_data["BF(dB)"])
			
			bayPass_data["gene"] = gea_data["gene"]
			bayPass_data["maf"] = gea_data["maf"]
			bayPass_data["pbar_qbar"] = gea_data["pbar_qbar"]

## Make a dataframe with all the PVEs per gene
		LA = gea_data.groupby(["gene"])["LA"].mean().to_frame()
## Make a dataframe with all the positions per gene (average position across SNPs)
		position = gea_data.groupby(["gene"])["position"].mean().to_frame()
## Make a dataframe with all the distances to the nearest PVE gene
		distance = distance_to_PVE(LA.to_dict()["LA"], position.to_dict()["position"], PVE_threshold)

## Do the top-candidate and WZA on uncorrected data AND do WZA on BayPass
		if bayPass != "none":
			wza_results_kendall = WZA( gea_data ,  "geno_k_tau_p_value", varName = "Z_kendall") 
			wza_results_bayPass = WZA( bayPass_data , "pVal", varName = "Z_bayPass") 
			TC_results_baypass = top_candidate( bayPass_data, "pVal" , label = "BP") 
			TC_results_kendall = top_candidate( gea_data, "geno_k_tau_p_value") 
 
			SNP_results_bayPass = simple_classifier( bayPass_data , "pVal" , bayPass = True, label = "BF" )
			SNP_results_kendall = simple_classifier( gea_data , "geno_k_tau_p_value"  ,label = "kendall")
## Do the top-candidate and WZA on uncorrected data
		else:
			TC_results = top_candidate( gea_data, "geno_k_tau_p_value") 
			wza_results = WZA( gea_data ,  "geno_k_tau_p_value") 
			SNP_results = simple_classifier( gea_data , "geno_k_tau_p_value" )


		if bayPass != "none":
			results = pd.concat([ distance, wza_results_kendall, wza_results_bayPass, TC_results_baypass, TC_results_kendall, SNP_results_bayPass, SNP_results_kendall, LA, position ] , axis = 1, sort = True ).reset_index()
			print(results)
		else:
			results = pd.concat([ distance, wza_results, TC_results, SNP_results, LA, position ] , axis = 1, sort = True ).reset_index()

		results["rep"] = rep

		yield( results )

""" R Code for the WZA would look something like:
gea$z <- qnorm(gea$pop_k_tau_p_value, lower.tail = F)
gea$z[gea$z==-Inf] <- qnorm(0.999, lower.tail = F)

gea_filt<-gea#[gea$maf > 0.01, ]

gene_pos <- tapply( gea_filt$position, as.factor(gea_filt$gene), mean)

LA <- tapply( gea_filt$LA, as.factor(gea_filt$gene), mean)

weiZ_num <- tapply( gea_filt$pbar_qbar * gea_filt$z, as.factor(gea_filt$gene), sum)

weiZ_den <- sqrt( tapply( gea_filt$pbar_qbar**2, gea_filt$gene, sum))

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
			type = str,
			dest = "bayPass", 
			help = "Do you want to include BayPass results? If so, give the directory here",
			default = "none")
	parser.add_argument("--directional", 
			action = "store_true",
			dest = "directional", 
			help = "Are the simulations of directional selection?")			
	parser.add_argument("--IBA", 
			action = "store_true",
			dest = "IBA", 
			help = "Do you want to restrict the analysis to the neutral fifth chromosome?",
			default = False)
	args = parser.parse_args()

## Make a list for each of the simulation parameter sets

	searchString = "/*.csv"

	simSet = [i for i in glob.glob(args.csv + searchString)]
		
	all_results = []

	for q in analyseSimSet(simSet, bayPass = args.bayPass, IBA = args.IBA):
		all_results.append( q.reset_index(drop=True) )
	myData = pd.concat( all_results, sort = True ).reset_index(drop=True)

	if "index" in list(myData):
		myData = myData.drop(["index"], axis =1)

	myData.to_csv( args.output, index = False)


if __name__ == "__main__":
	main()
