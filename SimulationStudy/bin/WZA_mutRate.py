import pandas as pd
import scipy.stats
import numpy as np
import sys, glob, argparse, random
from collections import Counter

from WZA import *



def WZA( GEA , statistic , MAF_filter = 0.05, varName = "", groupVar = "gene"):
## gea - the name of the pandas dataFrame with the gea results
## statistic - the name of the column with your p-values
## MAF_filter - the lowest MAF you wil tolerate
## varName - the column name for the weigted Z results

## Very small p-values throw Infinities when converted to z_scores, so I convert them to small numbers (i.e. 1e-15), or real big ones (1-1e-3)

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


def WZA_SNP_sampler( GEA , statistic , MAF_filter = 0.05, varName = "", groupVar = "gene"):
    num_SNP_list = np.array( [s[1].shape[0] for s in GEA.groupby(groupVar) ] )

    max_SNP_count = int(np.percentile(num_SNP_list[num_SNP_list!=0],50))

    WZA_output = []
    for gs in GEA.groupby(groupVar):
        if gs[1].shape[0] > max_SNP_count:
            TEMP = pd.concat( [WZA( gs[1].sample(max_SNP_count) ,  "empR_pVal", varName = "empR") for i in range(3)])
#            print(TEMP)
            new = pd.DataFrame(TEMP.mean()).transpose()
            new[groupVar] = gs[0]
            WZA_output.append( new.set_index(groupVar) )
    #        print (gs[0] )
    #        WZA_output.append(TEMP.mean(axis = 1))
    #        TEMP.loc[gs[0]] = TEMP.mean()
    #        print("hello", TEMP)
    #        WZA_output.append( TEMP.iloc[[-1]] )

        else:
            WZA_output.append( WZA( gs[1] ,  "empR_pVal", varName = "empR") )
#            print(WZA( gs[1] ,  "empR_pVal", varName = "empR"))

    return( pd.concat(WZA_output) )


def WZA_SNP_bins( GEA , statistic , MAF_filter = 0.05, varName = "", groupVar = "gene"):

    num_SNP_list = np.array( [s[1].shape[0] for s in GEA.groupby(groupVar) ] )

    intervals =  pd.IntervalIndex.from_breaks([0,10,100,1000])

    bin_names =np.array(["bin1","bin2","bin3"])

    wza_results_temp = WZA( GEA ,  "empR_pVal", varName = "empR")

    gene_bins =  [ bin_names[intervals.contains(n)][0] for n in num_SNP_list ]

    wza_results_temp["bin"] = gene_bins

    approx_pVals = []
    for gs in wza_results_temp.groupby("bin"):
        mean_Z = gs[1]["empR_Z"].mean()
        sd_Z = gs[1]["empR_Z"].std()
        temp =  -1*np.log10( 1-scipy.stats.norm.cdf( gs[1]["empR_Z"], mean_Z, sd_Z) )
        new = gs[1].copy()
        new["empR_Z"] = temp
        approx_pVals.append( new )

    return( pd.concat(approx_pVals))

def combineMutationRates(CSVs, mutation_map):
    output = []
    for j, k in zip(CSVs, mutation_map):
        output.append(j[j["gene"].isin(k)])
    new_df = pd.concat(output)
    new_df.sort_values("position")
    return(new_df)

## A master function to analyse a directory full of simulation results
def analyseSimSetMutRate(simulationSets, PVE_threshold = 0.005, MAF_filter = 0.05, SNP_sampler=False, SNP_bins = False):
    print(len(simulationSets[2]))
## Iterate over all GEA files in the directory
    for r in range(30):
        simSet = simulationSets[r]
        rep = simSet[0].split("/")[-1].split("_")[0]
        print("Analysing replicate:",rep)
## Read in the uncorrected GEA

        gea_data_1 = pd.read_csv( simSet[0] )
        gea_data_2 = pd.read_csv( simSet[1] )
        gea_data_3 = pd.read_csv( simSet[2] )

        all_genes = ["gene"+str(i) for i in range(1000)]
        random.shuffle(all_genes)
        mutMap = [all_genes[0:333],
                all_genes[333:666],
                all_genes[666:1000]]

        gea_data = combineMutationRates( [gea_data_1,gea_data_2,gea_data_3], mutMap )


## Due to a typo, some of the sims had pbar_q_bar in the header instead of pbar_qbar, so I just fix that here
        if "pbar_q_bar" in list(gea_data):
            gea_data.rename(columns={"pbar_q_bar":"pbar_qbar"}, inplace = True)


## Perform the MAF filter
        maf_filt_bool =  gea_data["maf"] > MAF_filter

        gea_data = gea_data[maf_filt_bool].copy()
## Calculate the empirical p-values on the correlation coefficient (Kendall's tau in this case)
        gea_data["empR_pVal"] =  scipy.stats.rankdata( abs(gea_data["geno_k_tau_p_value"]) )/ len(gea_data["geno_k_tau_p_value"])

## Do the analysis on Kendall's tau
        SNP_results = simple_classifier( gea_data , "geno_k_tau_p_value" ,label = "kendall")



## Make a dataframe with all the PVEs per gene
        LA = gea_data.groupby(["gene"])["LA"].mean().to_frame()
## Make a dataframe with all the positions per gene (average position across SNPs)
        position = gea_data.groupby(["gene"])["position"].mean().to_frame()
## Make a dataframe with all the distances to the nearest PVE gene
        distance = distance_to_PVE(LA.to_dict()["LA"], position.to_dict()["position"], PVE_threshold)

        if SNP_sampler:
            wza_results_corr = WZA_SNP_sampler( gea_data ,  "empR_pVal", varName = "empR")
#            print("!!!", wza_results_corr)
        elif SNP_bins:
            wza_results_corr = WZA_SNP_bins( gea_data ,  "empR_pVal", varName = "empR")
        else:
            wza_results_corr = WZA( gea_data ,  "empR_pVal", varName = "empR")

        wza_results_kendall = WZA( gea_data ,  "geno_k_tau_p_value", varName = "kendall")
        TC_results_kendall = top_candidate( gea_data, "geno_k_tau_p_value")

        results = pd.concat([ distance,
            wza_results_kendall,
            wza_results_corr,
            TC_results_kendall,
            SNP_results,
            LA,
            position ] , axis = 1, sort = True )
## Change this to make a sequential thing - like, do this by default then add BayPass then add RDA

        results["rep"] = rep

        results.reset_index(inplace = True)
        results.rename(columns = {"index":"gene"}, inplace = True)
        print(results)
        yield( results )

def main():
## Define command line args
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("--csv_m1",
        required = True,
        dest = "csv_m1",
        type = str,
        help = "The directory containing the CSVs for the first kind of mutation rate")
  parser.add_argument("--csv_m2",
        required = True,
        dest = "csv_m2",
        type = str,
        help = "The directory containing the CSVs for the second kind of mutation rate")
  parser.add_argument("--csv_m3",
        required = True,
        dest = "csv_m3",
        type = str,
        help = "The directory containing the CSVs for the third kind of mutation rate")
  parser.add_argument("--output",
        required = True,
        dest = "output",
        type = str,
        help = "The name of the output file")

  parser.add_argument("--sampleSNPs",
        action = "store_true",
        dest = "sampleSNPs",
        help = "Do you want to use the SNP sampling protocol?")

  parser.add_argument("--SNP_bins",
        action = "store_true",
        dest = "SNP_bins",
        help = "Do you want to use the SNP binning protocol?")

  args = parser.parse_args()

  if args.sampleSNPs == True and args.SNP_bins == True:
      return
## Make a list for each of the simulation parameter sets

  searchString = "/*.csv"

  simSet_m1 = [i for i in glob.glob(args.csv_m1 + searchString)]
  simSet_m2 = [i for i in glob.glob(args.csv_m2 + searchString)]
  simSet_m3 = [i for i in glob.glob(args.csv_m3 + searchString)]

  simSets = list(map(list, zip(*[ simSet_m1, simSet_m2, simSet_m3 ])))

  all_results = []

  for q in analyseSimSetMutRate(simSets,
                                SNP_sampler = args.sampleSNPs,
                                SNP_bins = args.SNP_bins):
      all_results.append( q.reset_index(drop=True) )

  myData = pd.concat( all_results, sort = True ).reset_index(drop=True)

  if "index" in list(myData):
      myData = myData.drop(["index"], axis =1)

  myData.to_csv( args.output, index = False)


if __name__ == "__main__":
  main()
