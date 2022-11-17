## A script to look at how the window size influences the WZA
from WZA import *
import argparse, glob
import pandas as pd
import numpy as np
## Read in the GEA results


def main():
## Define command line args
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--csv", 
			required = True,
			dest = "csv",
			type = str, 
			help = "The analysis file")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "The name of the output file")
	parser.add_argument("--nSNPs", "-n", 
			required = True, 
			dest = "nSNPs",
			type = int, 
			help = "The number of SNPs in an interval")
	args = parser.parse_args()
	
	n= args.nSNPs
	MAF_filter = 0.05


	searchString = "/*.csv"

	simSet = [i for i in glob.glob(args.csv + searchString)]
	gene_dataframes = []

	for gea_file in simSet:
		print(gea_file)
		rep = gea_file.split("/")[-1].split("_")[0]
		gea_data = pd.read_csv(gea_file)
		gea_data["empR_pVal"] =  1 - scipy.stats.rankdata( gea_data["geno_k_tau"]**2 )/ len(gea_data["geno_k_tau"])

		for gene in ["gene"+str(i) for i in range(1000)]:
	#	for gene in ["gene0"]:
		#	y = np.array([i for i in range(50)])
	# Make a copy of the results DF for each gene and reset the index
			gene_gea_data = gea_data[gea_data["gene"] == gene ].copy().reset_index(drop= True)
	## In this script, I have to apply the MAF filter early so that I can make the SNP sets
			gene_gea_data = gene_gea_data[gene_gea_data["maf"] > MAF_filter].copy().reset_index(drop= True)

			numSNPs = gene_gea_data.shape[0]
			print(gene)

	## Every n SNPs
			y = [l for l in range(numSNPs)]
			SNP_sets = np.array([[q,r] for q,r in zip(y[:numSNPs:n],y[n:numSNPs:n]) ])
			SNP_set_length_dict = {"set666_666":-99}
			for b in range(len(SNP_sets)):
				SNP_set_length_dict["set"+str(b)] = ( gene_gea_data["position"][SNP_sets[b][1]-1] - gene_gea_data["position"][SNP_sets[b][0]])

			# This horrible line makes a flattened list out of a list of lists for the indexes of the 
			SNP_set_labels = [ "set"+str(item) for sublist  in [n*[i] for i in range(SNP_sets.shape[0]) ] for  item in sublist]

	##  there may be orphan SNPs not assigned to any set. We'll call those set666_666 - those get ignored later on.
			SNP_set_labels += ["set666_666"]*(numSNPs - len( SNP_set_labels))
			gene_gea_data["SNP_set"] = SNP_set_labels

			gene_gea_data["SNP_set_interval"] = gene_gea_data["SNP_set"].map(SNP_set_length_dict)

	## Make a dataframe with all the PVEs per gene
			LA = gene_gea_data.groupby(["SNP_set"])["LA"].mean().to_frame()
	## Make a dataframe with all the positions per gene (average position across SNPs)
			position = gene_gea_data.groupby(["SNP_set"])["position"].mean().to_frame()
	## Make a dataframe with all the positions per gene (average position across SNPs)
			intervalWidths = gene_gea_data.groupby(["SNP_set"])["SNP_set_interval"].mean().to_frame()

#			TC_results = top_candidate( gea_data, "geno_k_tau_p_value") 
#			wza_results = WZA( gea_data ,  "geno_k_tau_p_value", varName = "kendall") 
#			wza_results_corr = WZA( gea_data ,  "empR_pVal", varName = "empR") 
#			SNP_results = simple_classifier( gea_data , "geno_k_tau_p_value" ,label = "SNP_kendall")

			wza_results_kendall = WZA( gene_gea_data ,  "geno_k_tau_p_value", varName = "kendall", groupVar = "SNP_set", MAF_filter = 0.05) 
			wza_results_empR = WZA( gene_gea_data ,  "empR_pVal", varName = "empR", groupVar = "SNP_set", MAF_filter = 0.05) 
			SNP_results_kendall = simple_classifier( gene_gea_data , "geno_k_tau_p_value"  ,label = "kendall", groupVar = "SNP_set")
			

			gene_results = pd.concat([ intervalWidths, wza_results_kendall, wza_results_empR,SNP_results_kendall, LA, position ] , axis = 1, sort = True ).rename(columns={"index":"SNP_set"})

			gene_results["gene"] = gene
			gene_results["rep"] = rep
#			print(genes_ranked)
			gene_dataframes.append(gene_results)

	final_output = pd.concat( gene_dataframes ).to_csv(args.output, index = False)
	
if __name__ == "__main__":
	main()
