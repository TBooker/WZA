import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse
from scipy.stats import norm


def WZA( gea , statistic , MAF_filter = 0.05, unweighted=False):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## NOTE, this function assumes that the DataFrame has a column named pbar_qbar
#	print("before sample\n", gea)

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	gea[statistic] = gea[statistic].clip(lower = 1e-15 )

	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

# convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z_score"] = scipy.stats.norm.ppf(1 - np.array( gea[statistic], dtype = float))

	gea["pbar_qbar"] = gea["MAF"]*(1-gea["MAF"])
	if unweighted == True:
		gea["pbar_qbar"] = 1

## Apply the MAF filter
	gea_filt = gea[ gea["MAF"] > MAF_filter ].copy()

	if gea_filt.shape[0]==0:
		return np.nan
## Calculate the numerator and the denominator for the WZA

	gea_filt["weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z_score"]

	gea_filt["weiZ_den"] =  gea_filt["pbar_qbar"]**2

	numerator = gea_filt["weiZ_num"].sum()

	denominator = np.sqrt(gea_filt["weiZ_den"].sum())

## We've calculated the num. and the den., let's take the ratio
	weiZ  = numerator/denominator

## Return the final dataframe
	return weiZ


def top_candidate( gea, thresh):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit - this should be the quantile threshold
##top_candidate_threshold - the probability point for calculating the expected number of genes
## Identifty the hits
	hits = ( gea["pVal"] < thresh ).sum()

	snps = gea.shape[0]

	top_candidate_p = scipy.stats.binom_test(hits, snps, thresh, alternative = "greater" )


	return top_candidate_p, hits

def adjust_WZA_with_spline(wza_df, roller = 50):

	# remove null Z values - they won't help us
	wza_t =  wza_df[ ~wza_df.Z.isnull() ].reset_index()
	wza_s = wza_t.sort_values('SNPs')

	rolled_Z_vars = wza_s.Z.rolling(window = roller).var()
	masking_array = ~rolled_Z_vars.isnull()
	rolled_Z_sd = np.sqrt(rolled_Z_vars)[masking_array]
	rolled_Z_means = wza_s.Z.rolling(window = roller).mean()[masking_array]
	rolled_mean_SNP_number = wza_s.SNPs.rolling(window = roller).mean()[masking_array]


	# Generating weights for polynomial function with degree =2 - standard deviation
	sd_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_sd, deg = 2)
	# Generating model with the given weights - the model as a function
	sd_polynomial_model = np.poly1d(sd_weights)

	# Generating weights for polynomial function with degree =2 - mean
	mean_weights = np.polyfit(rolled_mean_SNP_number, rolled_Z_means, deg = 2)
	# Generating model with the given weights - the model as a function
	mean_polynomial_model = np.poly1d(mean_weights)

#	xp = np.linspace(1,wza["SNPs"].max(),100)

	sd_predictions = sd_polynomial_model( wza_df["SNPs"] )
	mean_predictions = mean_polynomial_model( wza_df["SNPs"] )


	wza_p_values = [1- norm.cdf(wza_df["Z"][i], loc = mean_predictions[i], scale = sd_predictions[i]) for i in range(wza_df.shape[0])]

	# I think I have to loop as I don't think stats.cdf allows vectorised input
	# for i in range(wza.shape[0]):
	#     p_val = 1- norm.cdf(wza["Z"][i], loc = mean_predictions[i], scale = sd_predictions[i])
	#     print(i,wza["SNPs"][i], wza["Z"][i],  mean_predictions[i], sd_predictions[i], p_val)
	#     wza_p_values.append(p_val)

	wza_df["Z_pVal"] = wza_p_values

	return wza_df


def main():

## Define command line args
	parser = argparse.ArgumentParser(description="A script that implements the WZA, a method for combining evidence across closely linked SNPs in GEA studies.")

	parser.add_argument("--correlations", "-c",
			required = True,
			dest = "correlations",
			type = str,
			help = "The file containing the correlations")

	parser.add_argument("--summary_stat", "-s",
			required = True,
			dest = "summary_stat",
			type = str,
			help = "The name of the column you are analysing")

	parser.add_argument("--window", "-w",
			required = True,
			dest = "window",
			type = str,
			help = "The name of column containing the windows you want to analyse")

	parser.add_argument("--output",
			required = True,
			dest = "output",
			type = str,
			help = "The name of the output file")

	parser.add_argument("--sample_snps",
			required = False,
			dest = "sample_snps",
			type = int,
			help = "[OPTIONAL] Give the number of SNPs you want to downsample to. Give -1 if you want to use the 75th percentile of the number of SNPs. Note that calculating the median within the script is slow, so you may want to run a dummy analysis, get the median number of SNPs then use that explicitly.",
			default = 0)

	parser.add_argument("--snp_window",
			required = False,
			dest = "snp_window",
			type = int,
			help = "[OPTIONAL] Give the number of SNPs you want to downsample to. Give -1 if you want to use the 75th percentile of the number of SNPs. Note that calculating the median within the script is slow, so you may want to run a dummy analysis, get the median number of SNPs then use that explicitly.",
			default = 0)

	parser.add_argument("--resamples",
			required = False,
			dest = "resamples",
			type = int,
			help = "[OPTIONAL] Give the number of times you want to resample WZA scores when the number of SNPs exceeds the sample_snps threshold. [100]",
			default = 100)

	parser.add_argument("--min_snps",
			required = False,
			dest = "min_snps",
			type = int,
			help = "[OPTIONAL] Give the minimum number of SNPs you are willing to analyse per window. Default is 3",
			default = 3)

	parser.add_argument("--large_i_small_p",
			required = False,
			action = "store_true",
			help = "[OPTIONAL] Give this flag if extreme values of the summary stat you're using are large values.")

	parser.add_argument("--top_candidate_threshold",
			required = False,
			dest = "top_candidate_threshold",
			type = float,
			help = "[OPTIONAL] Give the percentile threshold you want to use for the top-candidate test (from Yeaman et al 2016 - Science)",
			default = 99)

	parser.add_argument("--verbose", "-v",
			required = False,
			action = "store_true",
			help = "[OPTIONAL] Give this flag if you want to run the script in verbose mode")

	parser.add_argument("--MAF",
			required = False,
			dest = "MAF",
			type = str,
			help = "[OPTIONAL] Give the name of the MAF colummn - if it's not explicitly given in the input dataframe")

	parser.add_argument("--sep",
			required = False,
			dest = "sep",
			type = str,
			default = "\t",
			help = "What separator do you use in your file?")

	parser.add_argument("--retain",
			required = False,
			dest = "retain",
			nargs="+",
			type = str,
			help = "Give any columns that you want to add to the output file")

	args = parser.parse_args()

	csv = pd.read_csv(args.correlations, sep = args.sep)
	if args.window not in list(csv):
		print("The window variable you provided is not in the dataframe you gave")
		return
#	if args.summary_stat not in list(csv):
#		print("The summary statistic variable you provided is not in the dataframe you gave")
#		return


	if "MAF" in list(csv):
		pass
	else:
		csv["MAF"] = csv[args.MAF].copy()

	if args.summary_stat =="parametric":
		csv["pVal"] = csv["geno_k_tau_p_value"]
	else:
		if args.large_i_small_p:
			if args.summary_stat =="RDA":
				csv["pVal"] =  1 - (csv[args.summary_stat]**2).rank()/csv.shape[0]
			else:
				csv["pVal"] =  1 - csv[args.summary_stat].rank()/csv.shape[0]
		else:
			csv["pVal"] = csv[args.summary_stat].rank()/csv.shape[0]


	if args.snp_window !=0:
		window_dat = []
		for win in set(	csv[args.window]):
			slice = csv[csv[args.window]==win].copy().reset_index()
			nrows = slice.shape
			ind = 0
			for s in range(0,slice.shape[0], args.snp_window):
				ind +=1
				if s+args.snp_window >slice.shape[0]:continue
				lil_slice = slice.iloc[s:s+args.snp_window].copy()
				lil_slice[args.window] =lil_slice[args.window]+"_"+str(ind)
				window_dat.append(lil_slice)
		csv = pd.concat(window_dat)


	csv_genes = csv[csv[args.window]!="none"]

	if args.verbose:
		print("here's a peek at the input data")
		print(csv_genes.head)


	csv_gb_gene = csv_genes.groupby(args.window)

	if args.sample_snps == -1:
## Get a list of the number of SNPs per contig

		csv_gb_gene_SNP_count = csv_genes.groupby(args.window)

		num_SNP_list = np.array( [s[1].shape[0] for s in csv_gb_gene_SNP_count if s[1].shape[0] > args.min_snps]  )

## Calculate the 75th percentile of SNPs per gene
		max_SNP_count = int(np.percentile(num_SNP_list[num_SNP_list!=0], 75))

		if args.verbose:
			print("Using the 75th percentile number of SNPs as the maximum in each gene:", max_SNP_count)

	elif args.sample_snps == 0:
		max_SNP_count = int(1e6) # This is just a large number that is never going to be the number of SNPs within a gene
		if args.verbose:
			print("The maximum number of SNPs in each gene:", max_SNP_count)

	else:
		max_SNP_count = int(args.sample_snps)
		if args.verbose:
			print("The maximum number of SNPs in each gene:", max_SNP_count)


	all_genes = []

	count = 0

	for g in csv_gb_gene:
		count += 1
		gene = g[0]
		gene_df = g[1].copy()

## Perform the WZA on the annotations in the contig using parametric ps

		if gene_df.shape[0] <=max_SNP_count:
## Perform the WZA on the annotations in the contig using empirical ps
			wza = WZA(gene_df, "pVal")
		else:
			wza = np.array( [ WZA(gene_df.sample(max_SNP_count), "pVal") for i in range(args.resamples ) ] ).mean()

#		print(wza)

		top_candidate_p, hits = top_candidate( gene_df, 1-(args.top_candidate_threshold/100))

#		print(top_candidate_p, hits, gene_df.shape[0], 1-(args.top_candidate_threshold/100))
#		print(wza_emp_pVal)

		if args.verbose:
			print("\ngene #:", count, "\tgene:", gene,"\tWZA:",wza,"\tTC:", top_candidate_p)

		output = {"gene": gene,
					"SNPs":gene_df.shape[0],
					"hits":hits,
					"gene":gene,
					"Z":wza,
					"top_candidate_p":top_candidate_p
					}

		all_genes.append( output )
#			input("\nPress Enter to continue...")



	if args.retain != None:
		WZA_DF_temp =  pd.DataFrame( all_genes )

		if args.verbose:
			print("\nAdding retained columns to the final dataframe")
		retained_df_list = []
		for r in args.retain:
			retained_df_list.append( csv.groupby(args.window)[r].mean() )
		retained_df = pd.concat(retained_df_list, axis = 1)
		WZA_DF_tmp =  pd.concat( [ WZA_DF_temp.set_index(args.window), retained_df ] , axis = 1).reset_index()

		WZA_DF = adjust_WZA_with_spline( WZA_DF_tmp )
#		output.to_csv(args.output, index = False)

	else:
		WZA_DF_tmp =  pd.DataFrame( all_genes )
		WZA_DF = adjust_WZA_with_spline( WZA_DF_tmp )

	WZA_DF.to_csv(args.output, index = False)






main()
