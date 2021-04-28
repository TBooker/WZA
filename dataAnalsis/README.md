## Analysis of data from Yeaman et al (2016)

I downloaded the data from Yeaman et al and stored it here ```../Yeaman_Data/```.

I grabbed out the big tables of individual SNP scores (Bayes factors or p-values), converted scores for DD_0 to empirical p-values and wrote them to this directory.

```
Rscript bin/getEmpiricalP_BF.R
```

<<<<<<< HEAD
The next step is to perform the WZA and the top candidate test to each contig in the dataframe:
=======

## NOTE: The above is deprecated as we are now no using the collapsed elements

For each set of GEA analyses performed on individual SNPs, we need to calculated the empirical pValue. The empirical p-value is just the rank of each SNP divided by the total number of SNPs. E.g. the smallest p-value would be 1/k, where k is the total number of SNPs tested. 

The following little script adds empirical p-values to the analysis files:

```
python bin/calculateEmpiricalP.py --correlations file_of_per_SNP_correlation_stats.txt --output file_of_per_SNP_correlation_stats_empiricalP.txt

```

The ```calculateEmpiricalP.py``` will make a dataframe for each environment, with empirical p-Values added as an additional column.

Once you've added the empirical-p values to the individual SNP files, you can run the WZA on them. Here's that step for analysing MAT:

```
python bin/WZA_CoAdapTree.py --correlations MAT_file_of_per_SNP_correlation_stats_empiricalP.txt --annotations DFtransv20202_DFrefedit.filtered.genes.collapsed.bed --output DF_interior_AHM.csv --env MAT --bed
```









>>>>>>> origin/master
