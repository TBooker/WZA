## Analyse GEA data, calculating Power and FDR
import pandas as pd
import argparse, glob
import numpy as np
import statsmodels.stats.multitest as multi
import scipy.stats as st

def grab_hits_SNPs(df, stat, threshold, large_i_small_p = False):
    if large_i_small_p == False:
        all_hits = (df[df[stat]<threshold]).gene
    elif large_i_small_p == True:
        all_hits = (df[df[stat]>threshold]).gene
    hit_genes = list(set(all_hits))
    true_positives = list(set((df[df["LA"]>0.005]).gene ))
    true_hits =  [p for p in hit_genes if p in true_positives]

    if len(hit_genes) == 0:
        power = 0
        fdr = np.nan
    else:
        power = len(true_hits)/len(true_positives)
        fdr = ( len(hit_genes) - len(true_hits) ) / len(hit_genes)

    return(power,fdr)


def main():
## Define command line args
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("--dir",
            required = True,
            dest = "dir",
            type = str,
            help = "The directory containing the files")
    parser.add_argument("--output",
            required = True,
            dest = "output",
            type = str,
            help = "The name for the output file")
    args = parser.parse_args()


    output = []
    for d in glob.glob("*output*/"):
        metadat  = d.split("/")[-2].split("_")
        map = metadat[0]
        demes = int(metadat[-2][1:])
        inds = int(metadat[-1][1:])
        if "mutMap" in metadat:
            mutation = "variation"
        else:
            mutation = "constant"
        for f in glob.glob(d+"/*singleSNPtests_wBP.csv"):
            rep = f.split("/")[-1].split("_")[0]


            csv = pd.read_csv(f)
            csv.dropna(axis = 0, inplace = True)
            # BayPass
            bp_res = grab_hits_SNPs(csv,
                "BF(dB)",
                20,
                large_i_small_p = True)
            output.append(["BF",map, mutation, rep, demes, inds, bp_res[0], bp_res[1]])

            lea_pass, lea_fdr = multi.fdrcorrection(csv["LEA"].values )
            csv["LEA_fdr"] = lea_fdr
            lea_res = grab_hits_SNPs(csv,
                "LEA_fdr",
                0.05,
                large_i_small_p = False)
            output.append(["LEA",map, mutation, rep, demes, inds, lea_res[0], lea_res[1]])

            mean_RDA = csv["RDA"].mean()
            sd_RDA = np.sqrt(csv["RDA"].var())

            csv["RDA_p"] =2*st.norm.cdf(-1.*abs(csv["RDA"]),
                                        loc = mean_RDA,
                                        scale = sd_RDA)
            rda_pass, rda_fdr = multi.fdrcorrection(csv["RDA_p"].values )

            csv["RDA_fdr"] = rda_fdr
            rda_res = grab_hits_SNPs(csv,
                "RDA_fdr",
                0.05,
                large_i_small_p = False)
            output.append(["RDA",map, mutation, rep, demes, inds, rda_res[0], rda_res[1]])

            wza = pd.read_csv(d+"/"+rep+"_WZA_geno_k_tau_p_value.csv")
            wza.dropna(axis = 0, inplace = True)

            wza_pass, wza_fdr = multi.fdrcorrection(wza["Z_pVal"].values )
            wza["WZA_fdr"] = wza_fdr

            wza_res = grab_hits_SNPs(wza,
                "WZA_fdr",
                0.05,
                large_i_small_p = False)
            output.append(["WZA",map, mutation, rep, demes, inds, wza_res[0], wza_res[1]])

    output_df = pd.DataFrame(output)
    output_df.columns = ["method","map","mutation","rep","demes","inds","power","fdr"]
    output_df.to_csv(args.output, index = False)
    print(output_df)

main()
