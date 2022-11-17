import argparse, pandas as pd
import numpy as np
from sklearn import metrics
import glob
### Here we'll write a script to summarise the unccorrected and WZA results into a single file for each map and set of conditions.
### The RScript I wrote was too janky and horrible to be a workable thing


def calcRecall(df_slice, numTrue, LA_threshold):
	if numTrue == 0:
		return np.nan
	elif numTrue > 0:
		genesFound = df_slice[df_slice["LA"] > LA_threshold ]
		return len(set(genesFound["gene"]))/numTrue

def calcPrecision(df_slice, numTrue, LA_threshold):
#	print("!",df_slice,numTrue)
	if numTrue == 0 or df_slice.shape[0] == 0 :
		return np.nan
	elif numTrue > 0:
		genesFound = df_slice[df_slice["LA"] > LA_threshold ]
		return len(set(genesFound["gene"]))/df_slice.shape[0]


def calcPVE(df_slice, cov_phenEnv):
	PCE_from_genesFound = sum(df_slice[df_slice["LA"] > 0 ]["LA"])
	return PCE_from_genesFound/cov_phenEnv


def calcFPR(df_slice, numFalse, LA_threshold):
	genesFound = df_slice[df_slice["LA"] > LA_threshold ]
	return( (df_slice.shape[0] - len(set(genesFound["gene"])))/numFalse )



def simple_classifier( gea, statistic, MAF_filter = 0.05, label = "", groupVar = "gene"):
	gea_filt = gea[ gea["maf"] > MAF_filter ].copy()
	## Get the maximum test statistic per gene
	## A little Pandas magic...
		# Group by "gene"
		# Get the maximum test statistic per gene
	if statistic in ["RDA","BF(dB)"]:
		if statistic == "RDA":
			gea_filt["RDA_sq"] = gea_filt["RDA"] **2
			gea_filt["empirical_p_raw"] = 1-(gea_filt["RDA_sq"].rank(method = "min")/ gea_filt.shape[0])
		else:
			gea_filt["empirical_p_raw"] = 1-(gea_filt[statistic].rank(method = "min")/ gea_filt.shape[0])
	else:
		gea_filt["empirical_p_raw"] = (gea_filt[statistic].rank(method = "min")/ gea_filt.shape[0])

# For each gene, just retain the SNP with the "most significant" hit
	genes_strimmed = gea_filt.loc[gea_filt.groupby('gene').empirical_p_raw.idxmin()]
	genes_strimmed["empirical_p"] = (genes_strimmed["empirical_p_raw"].rank(method = "min")/ genes_strimmed.shape[0])

	return genes_strimmed



def Precision_Recall(csv, LA_thresh, testStat = False, WZA_file = False):


	r = csv.split("/")[-1].split("_")[0]
	print("replicate",r)


## Calculate the empirical p-value for the test
	if WZA_file == True:
		temp = pd.read_csv(csv)

		Hits = temp["LA"] > LA_thresh

		hitGenes = set(temp[Hits]["gene"])

		numberOfHits = len(hitGenes)

		numFalse = temp.shape[0] - len(hitGenes)

		if testStat == False:
			temp["empirical_p"] = (temp["Z_pVal"].rank(method = "min")/ temp.shape[0])
			file_name = csv.split("/")[-1].strip(".csv")
			stat = "Z_"+file_name.replace(r+"_WZA_", "")
		elif testStat == "top_candidate_p":
			temp["empirical_p"] = (temp["top_candidate_p"].rank(method = "min")/ temp.shape[0])
			stat = "top_candidate_test"
	else:
		df = pd.read_csv(csv)
		if sum(df[testStat].isna()) == df.shape[0]:
			return None, None
		relevant_slice = df[["gene", "maf","LA", testStat]]

		temp = simple_classifier( relevant_slice, testStat, MAF_filter = 0.05, label = "slice", groupVar = "gene")

		Hits = temp["LA"] > LA_thresh

		stat = testStat

		hitGenes = set(temp[Hits]["gene"])

		numberOfHits = len(hitGenes)

		numFalse = temp.shape[0] - len(hitGenes)
		print( numberOfHits, numFalse)
	output = []
	for n in range(1,1001):
		if n%100==0:
			print("\t",n)
## Uncorrected data
		slice = temp[temp["empirical_p"] <= (n)/temp.shape[0]]

		output.append( {"rep":r,
					"top":n,

					"Recall":calcRecall(slice, numberOfHits, LA_thresh),

					"Precision":calcPrecision(slice, numberOfHits, LA_thresh),

					"FPR":calcFPR(slice, numFalse, LA_thresh)} )

	PR_data = pd.DataFrame(output)

	AUC_PR = metrics.auc(PR_data.Recall, PR_data.Precision)

	RESULT_dict = {"rep":r,
	"auc_pr":AUC_PR,
	"stat":stat}
	print(RESULT_dict)
	return( PR_data, RESULT_dict )

def summarizeResults(outputDF, bayPass_flag = False, rda_flag = False):
	out_dicts = []
	for r in set(outputDF["rep"]):
		slice = outputDF[outputDF["rep"]==r]

#		print(list(slice))

		out_DF = {"rep":r,
				"TC_uncor": metrics.auc(slice.Recall_TC_uncor, slice.Precision_TC_uncor),
				"SNP_uncor": metrics.auc(slice.Recall_SNP_uncor, slice.Precision_SNP_uncor),
				"WZA_uncor": metrics.auc(slice.Recall_WZA_uncor, slice.Precision_WZA_uncor),
				"WZA_emp": metrics.auc(slice.Recall_WZA_empR, slice.Precision_WZA_empR)}
		if bayPass_flag:
			out_DF["TC_bayP"] = metrics.auc(slice.Recall_TC_bayP, slice.Precision_TC_bayP)
			out_DF["SNP_bayP"] = metrics.auc(slice.Recall_SNP_bayP, slice.Precision_SNP_bayP)
			out_DF["WZA_bayP"] = metrics.auc(slice.Recall_WZA_bayP, slice.Precision_WZA_bayP)
		if rda_flag:
			out_DF["TC_RDA"] = metrics.auc(slice.Recall_TC_RDA, slice.Precision_TC_RDA)
			out_DF["SNP_RDA"] = metrics.auc(slice.Recall_SNP_RDA, slice.Precision_SNP_RDA)
			out_DF["WZA_RDA"] = metrics.auc(slice.Recall_WZA_RDA, slice.Precision_WZA_RDA)

		out_dicts.append(out_DF)

	return out_dicts

	"""
summarise_WZA <- function(top_hits){
  wza <-  as.data.frame(aggregate( truePositive_WZA ~ top, data = top_hits,
                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)))) )
  wza[,2:3] = unlist(wza$truePositive_WZA)
  names( wza ) <- c("top", "mean", "se")
  wza$test <- "WZA"
  tc <-  aggregate( truePositive_TC ~ top, data = top_hits,
                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
  tc[,2:3] = unlist(tc$truePositive_TC)
  names( tc ) <- c("top", "mean", "se")
  tc$test <- "Top-Candidate"
  snp <- aggregate( truePositive_SNPs ~ top, data = top_hits,
                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
  snp[,2:3] = unlist(snp$truePositive_SNPs)
  names( snp ) <- c("top", "mean", "se")
  snp$test <- "SNP-based"

  return( rbind( wza, tc, snp) )
}
"""

def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--result",
		required = True,
		dest = "result",
		type = str,
		help = "The directory containing the wza results")

	parser.add_argument("--singleSNP",
		required = True,
		dest = "singleSNP",
		type = str,
		help = "The directory containing the single SNP results")

#	parser.add_argument("--summary_stat",
#		required = True,
#		dest = "summary_stat",
#		type = str,
#		help = "The column from the results you want to analyse")

	parser.add_argument("--LA",
		required = True,
		dest = "LA",
		type = float,
		help = "Give the threshold for calling something a hit")

	parser.add_argument("--baypass",
		action = "store_true",
		dest = "baypass",
		help = "Give this flag if the single SNP results contain BayPass")

	parser.add_argument("--output",
		required = True,
		dest = "output",
		type = str,
		help = "Give the name for the output file. two files will be produced, so don't give the file extension!")


	args = parser.parse_args()

	wza_files = sorted([g for g in glob.glob(args.result+"/*WZA*csv")])

	single_snp_stats = ["geno_k_tau_p_value", "LEA", "RDA"]

	if args.baypass:
		snp_files = sorted([g for g in glob.glob(args.singleSNP+"/*singleSNPtests_wBP.csv")])
		single_snp_stats.append("BF(dB)")
	else:
		snp_files = sorted([g for g in glob.glob(args.singleSNP+"/*singleSNPtests.csv")])
		print(snp_files)

	pr_datasets = []
	result_dicts = []

	for s in snp_files:

#		print(s)
		for snp_stat in single_snp_stats:
#			print("hello")
#			relevant_slice = csv[["gene", "maf","LA",snp_stat]]
#			slice_classified = simple_classifier( relevant_slice, snp_stat, MAF_filter = 0.05, label = "slice", groupVar = "gene")
			pr_data, result_dict  = Precision_Recall(s,
													args.LA,
													testStat = snp_stat)
			if pr_data is None:
				continue
			pr_datasets.append(pr_data)
			result_dicts.append(result_dict)



	for w in wza_files:

		print(w,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		r = w.split("/")[-1].split("_")[0]
		file_name = w.split("/")[-1].strip(".csv")
		stat = file_name.replace(r+"_WZA_", "")
		if stat =="geno_k_tau_p_value":
			pr_data_TC, result_dict_TC  = Precision_Recall(w,
											args.LA,
											WZA_file = True,
											testStat = "top_candidate_p")
			print(pr_data_TC, "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
			if pr_data_TC is None:
				continue
			pr_datasets.append(pr_data_TC)
			result_dicts.append(result_dict_TC)

		pr_data, result_dict  = Precision_Recall(w,
											args.LA,
											WZA_file = True)
		print(pr_data, "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

		if pr_data is None:
			continue
		pr_datasets.append(pr_data)
		result_dicts.append(result_dict)
#		if count == 20:break

	print( pr_data )

	auc_results = pd.DataFrame(result_dicts)
	auc_results.to_csv("AUC_PrecisionRecall_"+args.output, index = False)

	pr_results = pd.concat(pr_datasets)
	pr_results.to_csv("PrecisionRecallCurves_"+args.output, index = False)

	return


main()
