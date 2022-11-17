import argparse, pandas as pd
import numpy as np
from sklearn import metrics
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

def Precision_Recall(csv, LA_thresh, BayPass = False, RDA = False, SNP_intervals = False):
	output_dicts = []
	print(set(csv["rep"]))
	for r in set(csv["rep"]):
#		if int(r) != 17:continue
		print("replicate",r)
		temp1 = csv[csv["rep"] == r].copy().reset_index()

		if SNP_intervals:
			temp2 = temp1[temp1["SNP_set_interval"] != -99].copy().reset_index()
			temp = temp2.iloc[temp2.groupby(["gene"])["empR_Z"].idxmax()].copy()

		else:
			temp = temp1.copy()

		Hits = temp["LA"] > LA_thresh
		hitGenes = set(temp[Hits]["gene"])
		numberOfHits = len(hitGenes)

		numFalse = temp.shape[0] - len(hitGenes)

## Calculate the empirical p-value for each of the tests
		temp["Z_empirical_p_uncor"] = 1 - (temp["kendall_Z"].rank(method = "min")/ temp.shape[0])

		temp["empZ_empirical_p_uncor"] = 1 - (temp["empR_Z"].rank(method = "min")/ temp.shape[0])

		if SNP_intervals:
			pass
		else:
			temp["TC_empirical_p_uncor"] = 1 - ((-1*np.log10(temp["top_candidate_p_TC2"])).rank(method = "first")/ temp.shape[0])

		temp["SNP_empirical_p_uncor"] =  (temp["kendall_rank"].rank(method = "min")/ temp.shape[0])

		if BayPass:
			temp["Z_empirical_p_bayP"] = 1 - (temp["bayPass_Z"].rank(method = "min")/ temp.shape[0])
			temp["TC_empirical_p_bayP"] = 1 - ((-1*np.log10(temp["BP_top_candidate_p_TC2"])).rank(method = "first")/ temp.shape[0])
			temp["SNP_empirical_p_bayP"] =  (temp["BF_rank"].rank(method = "min")/ temp.shape[0])

		if RDA:
			temp["Z_empirical_p_RDA"] = 1 - (temp["RDA_Z"].rank(method = "min")/ temp.shape[0])
			temp["TC_empirical_p_RDA"] = 1 - ((-1*np.log10(temp["RDA_top_candidate_p_TC2"])).rank(method = "first")/ temp.shape[0])
			temp["SNP_empirical_p_RDA"] =  (temp["RDA_rank"].rank(method = "min")/ temp.shape[0])

## Get the total Covariance(phen, env) across all genes - remember that I have to exclude adjacent genes (i.e. marked with LA = -99)
		Cov_phenEnv = sum(temp1[temp1["LA"] >0]["LA"])

		for n in range(1,1001):
			if n%100==0:
				print("\t",n)
## Uncorrected data
			WZA_slice_uncor = temp[temp["Z_empirical_p_uncor"] <= (n)/temp.shape[0]]
			WZA_slice_empR = temp[temp["empZ_empirical_p_uncor"] <= (n)/temp.shape[0]]



			if SNP_intervals:
				pass
			else:
				TC_slice_uncor = temp[temp["TC_empirical_p_uncor"] <= (n)/temp.shape[0]]

			SNP_slice_uncor = temp[temp["SNP_empirical_p_uncor"] <= (n)/temp.shape[0]]

			if BayPass:
				WZA_slice_bayPass = temp[temp["Z_empirical_p_bayP"] <= (n)/temp.shape[0]]
				TC_slice_bayPass = temp[temp["TC_empirical_p_bayP"] <= (n)/temp.shape[0]]
				SNP_slice_bayPass = temp[temp["SNP_empirical_p_bayP"] <= (n)/temp.shape[0]]

			if RDA:
				WZA_slice_RDA = temp[temp["Z_empirical_p_RDA"] <= (n)/temp.shape[0]]
				TC_slice_RDA = temp[temp["TC_empirical_p_RDA"] <= (n)/temp.shape[0]]
				SNP_slice_RDA = temp[temp["SNP_empirical_p_RDA"] <= (n)/temp.shape[0]]
#			print( calcFPR(WZA_slice_uncor, n, LA_thresh) )
#			print( calcFPR(SNP_slice_bayPass, n, LA_thresh))

#			print(calcPVE(WZA_slice_empR, Cov_phenEnv))
#			print(calcPVE(TC_slice_uncor, Cov_phenEnv))
#			print("")
			if SNP_intervals:
				output = {"rep":r,
					"top":n,
					"TP_WZA_uncor":calcTruePositives(WZA_slice_uncor, numberOfHits, LA_thresh),
					"TP_WZA_empR":calcTruePositives(WZA_slice_empR, numberOfHits, LA_thresh),
					"TP_SNP_uncor":calcTruePositives(SNP_slice_uncor, numberOfHits, LA_thresh),

					"FPR_WZA_uncor":calcFPR(WZA_slice_uncor, numFalse, LA_thresh),
					"FPR_SNP_uncor":calcFPR(SNP_slice_uncor, numFalse, LA_thresh),
					"FPR_WZA_empR":calcFPR(WZA_slice_empR, numFalse, LA_thresh)
					}
			else:
					output = {"rep":r,
								"top":n,

								"Recall_TC_uncor":calcRecall(TC_slice_uncor, numberOfHits, LA_thresh),
								"Recall_SNP_uncor":calcRecall(SNP_slice_uncor, numberOfHits, LA_thresh),
								"Recall_WZA_uncor":calcRecall(WZA_slice_uncor, numberOfHits, LA_thresh),
								"Recall_WZA_empR":calcRecall(WZA_slice_empR, numberOfHits, LA_thresh),

								"Precision_TC_uncor":calcPrecision(TC_slice_uncor, numberOfHits, LA_thresh),
								"Precision_SNP_uncor":calcPrecision(SNP_slice_uncor, numberOfHits, LA_thresh),
								"Precision_WZA_uncor":calcPrecision(WZA_slice_uncor, numberOfHits, LA_thresh),
								"Precision_WZA_empR":calcPrecision(WZA_slice_empR, numberOfHits, LA_thresh),

								"FPR_TC_uncor":calcFPR(TC_slice_uncor, numFalse, LA_thresh),
								"FPR_SNP_uncor":calcFPR(SNP_slice_uncor, numFalse, LA_thresh),
								"FPR_WZA_uncor":calcFPR(WZA_slice_uncor, numFalse, LA_thresh),
								"FPR_WZA_empR":calcFPR(WZA_slice_empR, numFalse, LA_thresh)}


			if BayPass:
				output["Recall_TC_bayP"] = calcRecall(TC_slice_bayPass, numberOfHits, LA_thresh)
				output["Recall_SNP_bayP"] = calcRecall(SNP_slice_bayPass, numberOfHits, LA_thresh)
				output["Recall_WZA_bayP"] = calcRecall(WZA_slice_bayPass, numberOfHits, LA_thresh)

				output["Precision_TC_bayP"] = calcPrecision(TC_slice_bayPass, numberOfHits, LA_thresh)
				output["Precision_SNP_bayP"] = calcPrecision(SNP_slice_bayPass, numberOfHits, LA_thresh)
				output["Precision_WZA_bayP"] = calcPrecision(WZA_slice_bayPass, numberOfHits, LA_thresh)

				output["FPR_TC_bayP"] = calcFPR(TC_slice_bayPass, numFalse, LA_thresh)
				output["FPR_SNP_bayP"] = calcFPR(SNP_slice_bayPass, numFalse, LA_thresh)
				output["FPR_WZA_bayP"] = calcFPR(WZA_slice_bayPass, numFalse, LA_thresh)

			if RDA:
				output["Recall_TC_RDA"] = calcRecall(TC_slice_RDA, numberOfHits, LA_thresh)
				output["Recall_SNP_RDA"] = calcRecall(SNP_slice_RDA, numberOfHits, LA_thresh)
				output["Recall_WZA_RDA"] = calcRecall(WZA_slice_RDA, numberOfHits, LA_thresh)

				output["Precision_TC_RDA"] = calcPrecision(TC_slice_RDA, numberOfHits, LA_thresh)
				output["Precision_SNP_RDA"] = calcPrecision(SNP_slice_RDA, numberOfHits, LA_thresh)
				output["Precision_WZA_RDA"] = calcPrecision(WZA_slice_RDA, numberOfHits, LA_thresh)

				output["FPR_TC_RDA"] = calcFPR(TC_slice_RDA, numFalse, LA_thresh)
				output["FPR_SNP_RDA"] = calcFPR(SNP_slice_RDA, numFalse, LA_thresh)
				output["FPR_WZA_RDA"] = calcFPR(WZA_slice_RDA, numFalse, LA_thresh)

			output_dicts.append( output )

	return( pd.DataFrame(output_dicts) )
#		length(droplevels(temp[temp$LA > LA_thresh,]$gene))

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

	parser.add_argument("--wza",
		required = True,
		dest = "wza",
		type = str,
		help = "The wza results that you want to summarise")

	parser.add_argument("--LA",
		required = True,
		dest = "LA",
		type = float,
		help = "Give the threshold for calling something a hit")

	parser.add_argument("--output",
		required = True,
		dest = "output",
		type = str,
		help = "Give the name for the output file. two files will be produced, so don't give the file extension!")

	parser.add_argument("--BayPass",
		dest = "BayPass",
		action = "store_true",
		help = "Give this flag if you are analysing BayPass results",
		default = False)

	parser.add_argument("--RDA",
		dest = "RDA",
		action = "store_true",
		help = "Give this flag if you are analysing RDA results",
		default = False)

	parser.add_argument("--SNP_intervals",
		dest = "SNP_intervals",
		action = "store_true",
		help = "Give this flag if you are analysing SNP interval results",
		default = False)

	args = parser.parse_args()

	wza = pd.read_csv(args.wza)

	outputSummary  = Precision_Recall(wza, args.LA, BayPass  = args.BayPass,  RDA  = args.RDA, SNP_intervals = args.SNP_intervals)
	outputSummary.to_csv(args.output + ".perRep.csv")
	AUC = summarizeResults(outputSummary, bayPass_flag = args.BayPass, rda_flag = args.RDA)
	pd.DataFrame(AUC).to_csv(args.output + ".AUC.csv")



main()
