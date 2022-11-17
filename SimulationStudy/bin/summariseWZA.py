import argparse, pandas as pd
import numpy as np
### Here we'll write a script to summarise the unccorrected and WZA results into a single file for each map and set of conditions. 
### The RScript I wrote was too janky and horrible to be a workable thing

def calcTruePositives(df_slice, numTrue, LA_threshold):
	if numTrue == 0:
		return np.nan
	elif numTrue > 0:
		genesFound = df_slice[df_slice["LA"] > LA_threshold ]
		return len(set(genesFound["gene"]))/numTrue
#	return sum( df_slice["LA"] > LA_threshold )

def calcPVE(df_slice, cov_phenEnv):
	PCE_from_genesFound = sum(df_slice[df_slice["LA"] > 0 ]["LA"])
	return PCE_from_genesFound/cov_phenEnv


def calcFDR(df_slice, nTotal, LA_threshold):
	genesFound = df_slice[df_slice["LA"] > LA_threshold ]
	return (nTotal - len(set(genesFound["gene"])))/nTotal

def NumberOfTruePositives(csv, LA_thresh, BayPass = False, SNP_intervals = False):
	output_dicts = []
	for r in set(csv["rep"]):
		print(r)
		temp1 = csv[csv["rep"] == r].copy()
		
		if SNP_intervals:
			temp2 = temp1[temp1["SNP_set_interval"] != -99].copy().reset_index()
			temp = temp2.iloc[temp2.groupby(["gene"])["empR_Z"].idxmax()].copy()

		else:
			temp = temp1.copy()

		Hits = temp["LA"] > LA_thresh
		hitGenes = set(temp[Hits]["gene"])

		numberOfHits = len(hitGenes)


## Calculate the empirical p-value for each of the tests
		temp["Z_empirical_p_uncor"] = 1 - (temp["kendall_Z"].rank()/ temp.shape[0])
		temp["empZ_empirical_p_uncor"] = 1 - (temp["empR_Z"].rank()/ temp.shape[0])
		if SNP_intervals:
			pass
		else:
			temp["TC_empirical_p_uncor"] = 1 - ((-1*np.log10(temp["top_candidate_p_TC2"])).rank()/ temp.shape[0])

		temp["SNP_empirical_p_uncor"] =  (temp["kendall_rank"].rank()/ temp.shape[0])
		if BayPass:
			temp["Z_empirical_p_bayP"] = 1 - (temp["bayPass_Z"].rank()/ temp.shape[0])
			temp["TC_empirical_p_bayP"] = 1 - ((-1*np.log10(temp["BP_top_candidate_p_TC2"])).rank()/ temp.shape[0])
			temp["SNP_empirical_p_bayP"] =  (temp["BF_rank"].rank()/ temp.shape[0])

## Get the total Covariance(phen, env) across all genes - remember that I have to exclude adjacent genes (i.e. marked with LA = -99)
		Cov_phenEnv = sum(temp1[temp1["LA"] >0]["LA"])

		for n in range(1,51):
			print(n)
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

#			print( WZA_slice_uncor )
#			print( calcFDR(WZA_slice_uncor, n, LA_thresh) )
#			print( calcFDR(SNP_slice_bayPass, n, LA_thresh))
			
#			print(calcPVE(WZA_slice_empR, Cov_phenEnv))
#			print(calcPVE(TC_slice_uncor, Cov_phenEnv))
#			print("")
			
			if BayPass:
				output = {"rep":r,
						"top":n,
						
						"TP_TC_uncor":calcTruePositives(TC_slice_uncor, numberOfHits, LA_thresh),
						"TP_TC_bayP":calcTruePositives(TC_slice_bayPass, numberOfHits, LA_thresh),
						"TP_SNP_uncor":calcTruePositives(SNP_slice_uncor, numberOfHits, LA_thresh),
						"TP_SNP_bayP":calcTruePositives(SNP_slice_bayPass, numberOfHits, LA_thresh),
						"TP_WZA_uncor":calcTruePositives(WZA_slice_uncor, numberOfHits, LA_thresh),
						"TP_WZA_bayP":calcTruePositives(WZA_slice_bayPass, numberOfHits, LA_thresh),
						"TP_WZA_empR":calcTruePositives(WZA_slice_empR, numberOfHits, LA_thresh),

						"PCE_TC_uncor":calcPVE(TC_slice_uncor, Cov_phenEnv),
						"PCE_TC_bayP":calcPVE(TC_slice_bayPass, Cov_phenEnv),
						"PCE_SNP_uncor":calcPVE(SNP_slice_uncor, Cov_phenEnv),
						"PCE_SNP_bayP":calcPVE(SNP_slice_bayPass, Cov_phenEnv),
						"PCE_WZA_uncor":calcPVE(WZA_slice_uncor, Cov_phenEnv),
						"PCE_WZA_bayP":calcPVE(WZA_slice_bayPass, Cov_phenEnv),
						"PCE_WZA_empR":calcPVE(WZA_slice_empR, Cov_phenEnv),


						"FD_TC_uncor":calcFDR(TC_slice_uncor, n, LA_thresh),
						"FD_TC_bayP":calcFDR(TC_slice_bayPass, n, LA_thresh),
						"FD_SNP_uncor":calcFDR(SNP_slice_uncor, n, LA_thresh),
						"FD_SNP_bayP":calcFDR(SNP_slice_bayPass, n, LA_thresh),
						"FD_WZA_uncor":calcFDR(WZA_slice_uncor, n, LA_thresh),
						"FD_WZA_bayP":calcFDR(WZA_slice_bayPass, n, LA_thresh),
						"FD_WZA_empR":calcFDR(WZA_slice_empR, n, LA_thresh)}
			else:
				if SNP_intervals:
					output = {"rep":r,
						"top":n,
						"TP_WZA_uncor":calcTruePositives(WZA_slice_uncor, numberOfHits, LA_thresh),
						"TP_WZA_empR":calcTruePositives(WZA_slice_empR, numberOfHits, LA_thresh),
						"TP_SNP_uncor":calcTruePositives(SNP_slice_uncor, numberOfHits, LA_thresh),
			
						"FD_WZA_uncor":calcFDR(WZA_slice_uncor, n, LA_thresh),
						"FD_SNP_uncor":calcFDR(SNP_slice_uncor, n, LA_thresh),
						"FD_WZA_empR":calcFDR(WZA_slice_empR, n, LA_thresh)
						}
				else:
					output = {"rep":r,
						"top":n,
						"TP_WZA_uncor":calcTruePositives(WZA_slice_uncor, numberOfHits, LA_thresh),
						"TP_WZA_empR":calcTruePositives(WZA_slice_empR, numberOfHits, LA_thresh),
						"TP_TC_uncor":calcTruePositives(TC_slice_uncor, numberOfHits, LA_thresh),
						"TP_SNP_uncor":calcTruePositives(SNP_slice_uncor, numberOfHits, LA_thresh),
			
						"FD_WZA_uncor":calcFDR(WZA_slice_uncor, n, LA_thresh),
						"FD_TC_uncor":calcFDR(TC_slice_uncor, n, LA_thresh),
						"FD_SNP_uncor":calcFDR(SNP_slice_uncor, n, LA_thresh),
						"FD_WZA_empR":calcFDR(WZA_slice_empR, n, LA_thresh)
						}

			output_dicts.append( output )
		
	return( pd.DataFrame(output_dicts) ) 
#		length(droplevels(temp[temp$LA > LA_thresh,]$gene))

def summarizeResults(outputDF):
	mean = outputDF.groupby(["top"]).mean()
	mean["rep"] = "mean"

	nReps = len(set(outputDF.rep)) 
	stdErr = outputDF.groupby(["top"]).std()/np.sqrt(nReps) 
	stdErr["rep"] = "stdErr"

	summaryDF = pd.concat([mean, stdErr]) 
	summaryDF.reset_index(level = 0, inplace = True)	
	return summaryDF
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

	parser.add_argument("--SNP_intervals", 
		dest = "SNP_intervals",
		action = "store_true",
		help = "Give this flag if you are analysing SNP interval results",
		default = False)

	args = parser.parse_args()

	wza = pd.read_csv(args.wza)
	
#	proportion_of_CVE(wza, BayPass  = args.BayPass, SNP_intervals = args.SNP_intervals)

#	return
	outputSummary  = NumberOfTruePositives(wza, args.LA, BayPass  = args.BayPass, SNP_intervals = args.SNP_intervals)
	outputSummary.to_csv(args.output + ".perRep.csv")
	summarizeResults(outputSummary).to_csv(args.output + ".summary.csv")

	


main()
