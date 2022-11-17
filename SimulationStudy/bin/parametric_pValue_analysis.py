import pandas as pd
import numpy as np
import scipy.stats as st

def significantSlice(DF, stat, threshold):
    sigSlice = DF[DF[stat]<threshold]
    genes = set(sigSlice.gene)
    LA = sigSlice.groupby(["gene"])["LA"].mean()
    hits = (LA > 0.005).sum()
    return hits, len(LA) - hits

def main():
    output = []
    for samp in [40,20,10]:
        if samp == 40:
            sample = ""
        elif samp == 20:
            sample = "_n20"
        elif samp == 10:
            sample = "_n10"
        for m in ["BC_Map", "cline", "trunc"]:
            if m == "BC_Map":
                mn = "1"
            elif m == "trunc":
                mn = "2"
            elif m == "cline":
                mn = "3"
            for s in ["0.0136", "0.003"]:
                wza = pd.read_csv( "~/work/GEA/simulations/SimulationStudy/s"+s+"/"+m+sample+".wza.csv" )


                for i in range(30):
                    rep = i+1
                    print("rep:", rep)
                    current = wza[wza["rep"] == rep].copy()


                    all_LA_genes = (current.groupby(["gene"])["LA"].mean() > 0.005).sum()

                    current["Z_pVal"] =  1- st.norm.cdf(current["empR_Z"],
                                                current["empR_Z"].mean(),
                                                np.sqrt(current["empR_Z"].var()))
                    gea = pd.read_csv("~/work/GEA/simulations/SimulationStudy/s"+s+"/"+m+"_output"+sample+"/" + str(rep) + "_"+s+"_"+mn+".12Loci.directionalSelection_d"+str(samp)+"n50.csv")

    #                rda = pd.read_csv("~/work/GEA/simulations/SimulationStudy/s"+s+"/"+m+sample+"_RDA/" + str(rep) + "_"+s+"_"+mn+".12Loci.directionalSelection_d40n50.rda.txt")

    #                rda["pVal"] =  (2*st.norm.cdf(abs(rda["RDA1"]),
    #                                            rda["RDA1"].mean(),
    #                                            np.sqrt(rda["RDA1"].var())))
    #                gea["rda_p"] = rda.pVal

            #        print( gea )

                    SNP_threshold = 1-(1-0.05)**(1/gea.shape[0])

                    window_threshold = 1-(1-0.05)**(1/1000)

                    empZ_tps, empZ_fps  = significantSlice(current, "Z_pVal",window_threshold )

                    TC_tps, TC_fps  = significantSlice(current, "top_candidate_p_TC2",window_threshold )

                    kendall_tps, kendall_fps  = significantSlice(gea, "geno_k_tau_p_value",SNP_threshold )

    #                rda_tps, rda_fps  = significantSlice(gea, "rda_p",SNP_threshold )

                    output.append( [s, m, samp,"WZA", empZ_tps/all_LA_genes, empZ_fps/(empZ_tps+ empZ_fps), rep ])
                    output.append( [s, m, samp,"TC", TC_tps/all_LA_genes, TC_fps/(TC_tps + TC_fps), rep ])
                    output.append( [s, m, samp,"kendall", kendall_tps/all_LA_genes, kendall_fps/(kendall_tps+ kendall_fps), rep ])
            #        print("\tRDA", rda_tps/all_LA_genes, rda_fps/(rda_tps + rda_fps) )


                print()
#        break
    output_DF = pd.DataFrame(output,
            columns = ["s","map","sampleSize","Test","Power","FDR","rep"])
    output_DF.to_csv("TEMPTMEPTEMP.csv", index = False)

    print(output_DF)
"""
        rda$pVal =  2*pnorm(abs(rda$RDA1),
                               mean(rda$RDA1),
                               sd(rda$RDA1),
                         lower.tail = F)
"""
main()
