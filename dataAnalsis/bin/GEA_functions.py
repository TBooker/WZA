import pandas as pd
import scipy.stats
import numpy as np
from GEA_functions import *

def WZA( gea , statistic , MAF_filter = 0.05, varName = "Z"):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results
## NOTE, this function assumes that the DataFrame has a column named pbar_qbar

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15)

	gea[statistic] = gea[statistic].clip(lower = 1e-15 )

	gea[statistic] = gea[statistic].replace(1 , 1-1e-3)

# convert the p-values into 1-sided Z scores (hence the 1 - p-values)
	gea["z"] = scipy.stats.norm.ppf(1 - np.array( gea[statistic], dtype = float))
## Apply the MAF filter
#	gea_filt = gea[ gea["maf"] > MAF_filter ].copy()

	gea_filt = gea.copy()

## Calculate the numerator and the denominator for the WZA

	gea_filt[varName+"_weiZ_num"] = gea_filt["pbar_qbar"] * gea_filt["z"]

	gea_filt[varName+"_weiZ_den"] =  gea_filt["pbar_qbar"]**2 	

	numerator = gea_filt.groupby(["gene"])[varName+"_weiZ_num"].sum().to_frame()

	denominator = np.sqrt(gea_filt.groupby(["gene"])[varName+"_weiZ_den"].sum()).to_frame()

## We've calculated the num. and the den., let's make a dataframe that has both 
	weiZ  = pd.concat([numerator,denominator], axis = 1, sort = False)

## Actually calculate the Z scores for each gene
	weiZ[varName] = weiZ[varName+"_weiZ_num"] / weiZ[varName+"_weiZ_den"]

## Deprecated
## One might be interested in calculating a p_value from the Z-scores (though this only works if the data are normal, which they won't be if there's populations structure).
##	weiZ[varName+"_hits"] = (weiZ[varName] > scipy.stats.norm.ppf(1 - 0.05/50)).astype(int)

## Return the final dataframe
	return weiZ


## A function for performing the top-candidate test

def top_candidate( gea, thresh, threshQuant, statistic, top_candidate_threshold, MAF_filter = 0.05):

## gea - the name of the pandas dataFrame with the gea results
## thresh - the p_value threshold for determining hits
## MAF_filter - the lowest MAF you wil tolerate
## prop_hits - the average probility of getting a hit - this should be the quantile threshold
##top_candidate_threshold - the probability point for calculating the expected number of genes
## Identifty the hits
	gea["hits"] = ( -np.log10( np.array( gea[statistic] , dtype = float)) > -np.log10(thresh)).astype(int)

## Apply the MAF filter
#	gea_filt = gea[ gea["maf"] > MAF_filter ]
	gea_filt = gea.copy()

## Count the hits per gene
	num_hits = gea_filt.groupby(["gene"])["hits"].sum().to_frame()

## Count the SNPs per gene
	num_SNPs = gea_filt.groupby(["gene"])["hits"].count().to_frame()

## Make a single DF with the hits and the SNPs
	TC  = pd.concat([num_hits, num_SNPs], axis = 1, sort = False) 

## Name the cols
	TC.columns = ["hits", "SNPs"]

	return(TC)

## Init an empty vector for p_values (the top-candidate index)
	p_vals = []

### Init an empty vector for expected hits at the "top-candidate" threshold
#	expectedHits = []
#	for index, row in TC.iterrows():
#		print(row.hits, row.SNPs, top_candidate_threshold )
#		p_vals.append( scipy.stats.binom_test(row.hits, row.SNPs, threshQuant, alternative = "greater" ) )
#		expectedHits.append( scipy.stats.binom.ppf( top_candidate_threshold , row.SNPs, threshQuant ) )


## Add the TC p_vals to the genes DF
#	TC["top_candidate_p"] = p_vals

## Add the TC expected number of hits to the genes DF
#	TC["expected_outliers"] = expectedHits


## Return the resulting dataFrame
#	return(TC)
