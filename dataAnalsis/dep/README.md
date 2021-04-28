Here's a little bit of test data to analyse CoAdapTree GEA results.

For each annotation in a specified file (either GFF or BED), the script ```bin/WZA_CoAdapTree.py``` performs both the WZA and the Top-Candidate test (Yeaman et al 2016) on the results from GEA.

To deal with overlapping annotations, I used BedTools merge as follows:
```
 ~/software/bedtools2/bin/bedtools merge -i DFtransv20202_DFrefedit.filtered.genes.gff -c 9 -o collapse > DFtransv20202_DFrefedit.filtered.genes.collapsed.bed
```
This takes a GFF file (in this case ```DFtransv20202_DFrefedit.filtered.genes.gff```) and merges overlapping elements. You could merge neighbouring elements using the ```-d``` command.

The script can be used to analyse a particular environment (from MAT MWMT MCMT TD MAP MSP AHM SHM DD_0 DDS NFFD bFFP FFP PAS EMT EXT Eref CMD), for example AHM can be analysed as follows:

```
python bin/WZA_CoAdapTree.py --correlations split_snp_env_spearmans_rho_interior_AHM_maf_TOM.txt --annotations DFtransv20202_DFrefedit.filtered.genes.collapsed.bed --output DF_interior_AHM.csv --env AHM --bed
```


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









