# WZA
A repository containing the analysis scripts and simulation configuration files to accompany our project on the Weighted-Z Analysis for genotype-environment association (GEA) studies.

We call the method the WZA for short, and it's pronounced the same way that Wu-Tang Clan fans might pronounce RZA or GZA or SZA.

A general purpose script for performing the WZA on input files (described below) can be found at [general_WZA_script.py](general_WZA_script.py). Instructions for usage of this script, along with a worked example are given below.

We also include a thoroughly worked through example here: [simulationDataVignette/](simulationDataVignette/)

This repository also contains code to simulate meta-populations adapting to heterogeneous environments and to perform GEA analyses on the resulting data. The directory [SimulationStudy/](SimulationStudy) contains more information on the simulations.  

If you are interested in performing the WZA for your own projects, feel free to get in touch if you run into any implementation issues - either by raising an issue on this repository, or via email (my email address is given in my GitHub profile).

____________________

## Performing the WZA on GEA results

We have implemented the WZA analysis as a Python (3.10) script [general_WZA_script.py](general_WZA_script.py). This script uses a number of third-party packages as well as some standard packages distributed with Python. Specifically, you will need the following packages installed and available to Python:
* pandas
* scipy
* numpy

These three packages are all standard Python packages that can easily be installed via ```pip```, ```conda``` or built form source.

The method was developed on Ubuntu and MacOS operating systems, and have been tested on both. If you are having difficulty running the script, feel free to get in touch (see above for contact info).

The script can be invoked at the command line easily:

```{sh}
python3 general_WZA_script.py

# Or, to see the help screen:
python3 general_WZA_script.py --help

```
________________________________________
#### Input GEA results

To perform the WZA as we did in the paper, you need to have first performed Kendall's tau-b (a rank correlation analysis that can handle tied datapoints) on each SNP that you are considering. Other summary statistics (such as squared loadings from RDA or Bayes factors from BayPass can easily be incorporated).

As input, the [general_WZA_script.py](general_WZA_script.py) requires one input file, a table of GEA scores for your SNPs. The script refers to this file as the "correlations" file.

This table is required to contain the following information:
* Scores for individual SNPs (e.g. these could be *p*-values from a correlation analysis or Bayes factors)
* A variable that indicates the windows that SNPs belong to (e.g. a column of gene names, or analysis window names)
* The average minor allele frequency across populations

With an input file (with header) that contains these things called "my_correlations.csv", the WZA can be performed simply using:


```{sh}
python3 general_WZA_script.py --correlations my_correlations.csv \
                              --summary_stat pVal \
                              --window window_names \
                              --MAF maf \
                              --output my_WZA.csv \
                              --sep ","
```
Let's break down what the above command actually does.

It runs the WZA on the file called "my_correlations.csv". It uses the values in the column called "pVal" as the basis of the WZA. It groups SNPs together on the basis of the column names "window_names". It uses the "maf" column as weights in the WZA analysis. It stores the output as "my_WZA.csv". Finally, because the input file was comma separated (as indicated by the file extension), we specified that as "," using the ```--sep``` argument.

Here are all the arguments that you can specify to the WZA script:

* ```--correlations, -c```: **The input file name**
* ```--summary_stat, -s```: **The statistic within your file that you want to use as the basis of the WZA**
* ```--window, -w```: **The column within the input file that indicates the windows to which SNPs belong**
* ```--output```: **The name of the output file you want to save**
* ```--sample_snps```: The percentile of the SNP number distribution that you want to set as the upper limit* [75]
* ```--resamples```: The number of resamples you want to perform for windows that have more SNPs than the upper limit set by ```sample_snps```[100]
* ```--min_snps```: The minimum number of SNPs you want to analyse for each window [3]
* ```--large_i_small_p```: A flag that indicates whether the ```summary_stat``` you are analyzing is descending or ascending in terms of evidence for GEA (e.g. *p*-values are descending, Bayes factors are ascending). [FALSE]
* ```--top_candidate_threshold```: The percentile for calling SNPs outliers from the genome-wide distribution (DOES NOT AFFECT THE WZA)[99]
* ```--verbose, -v```: Give this flag to specify you want to get output from the program as it runs
* ```--MAF```: The column containing the average minor allele frequency across populations
* ```--sep```: The field you use to separate columns in the input file ["\t"]
* ```--retain```: List any columns you want to retain from the input file in the output. Will return the mean of numeric values.

*Required arguments are given in bold*

*Default values for optional parameters are given in square brackets*

## Worked example

In our paper describing the WZA, we performed an analysis of some previously published data for the lodgepole pine (generated by Yeaman et al. 2016 - Science). Below we include a walkthrough of how we performed the WZA on that data using the [general_WZA_script.py](general_WZA_script.py).

The data that Yeaman et al. (2016) included with their paper contains GEA results for genome-wide SNPs for populations of lodgepole pine collected from across North Western North America. For each SNP, Yeaman et al. (2016) performed Spearman's rank correlation (i.e. Spearman's rho) for each SNP with numerous environmental variables. In the paper we only re-analysed one of these, DD0 - degrees days below 0°C (i.e. chilling degree-days).

In their paper, Yeaman et al. (2016) used a method they developed called the "top-candidate method". This method is essentially a binomial test that asks whether a particular genomic region has more outliers than expected given the genome-wide average. The [general_WZA_script.py](general_WZA_script.py) implements this method along with the WZA.
______________________

We downloaded the data from Yeaman et al. (2016) from [the DataDryad repository containing their data](https://datadryad.org/stash/dataset/doi:10.5061/dryad.0t407). We don't include the full file here as it is too big to be stored on GitHub. We do, however, include a file containing the first 10,000 SNPs from the file. That file is at this location: [testData/Yeaman_demoData.csv](testData/Yeaman_demoData.csv).

Here's how we analyse the data using the [general_WZA_script.py](general_WZA_script.py):

``` {sh}
# First grab a chunk of Yeaman et al's data...
head -n10001  var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_raw | gzip > Yeaman_demoData.tsv.gz

# That file that is in the Yeaman et al 2016 Dryad repo
# I gzipped the file to make GitHub easier
```

The file is really too big too walk through (there are 72 columns!). The file contains data for each SNP in their data set (we sampled the first 10,000 from these). For each SNP, Yeaman et al performed Spearman's rho for allele frueqency and numerous environmental variables. The p-value from those correlation tests are stored in the file under headers corresponding to each environmental variable.

So to perform the WZA on this data, we do the following:

``` {sh}


python3 general_WZA_script.py --correlations testData/Yeaman_demoData.tsv.gz \
                              --summary_stat DD_0 \
                              --window tcontig \
                              --MAF minor_freq \
                              --output testData/my_WZA.csv \
                              --sep "\t"

```

Running this command will produce a file that contains entries for 302 genes - Yeaman et al (2016) used contigs from the pine transcriptome as the windows in their analyses.

## Choosing window sizes

The WZA is an analysis window based genome scan. As with other analyses of this type, setting the width of analysis windows is something that requires careful consideration. In the Yeaman et al. (2016) analysis, the authors chose to use gene annotations as the windows when applying their top-candidate method, which we followed on from when re-analysing their data with the WZA.

Using gene annotations as windows in genome scans was done because the data the authors were working with was exome-capture data, so all SNPs were within or closely linked to genes. Using genes as analysis windows has a real benefit, the results for a particular window can be tied to a particular biological function pretty unambiguously. However, there are a couple of downsides to this approach. First, gene lengths can vary a great deal (over a couple of orders of magnitudes in conifers, for example). Such variation in window width has a similar effect on the distribution of genome-scan summary statistics as varying the recombination rate across windows. Specifically, longer windows potentially include more distinct evolutionary histories than shorter windows. This difference influences the variance in genome scan summary statistics, potentially leading to heteroskedasticity in genome scan results. Additionally, for other forms of sequence data, there may be many SNPs in intergenic regions, so using genes as windows effectively throws all of that data away.

In general, we recommend keeping the extend of LD in the genome in mind when setting the width of windows for the WZA. In our simulations validating the method, LD broke down to background levels within 2.5kbp or so for neutrally evolving genes. For genes that contained locally adaptive SNPs, LD broke down over longer distances. We used analysis windows of 10kbp in our analyses and found that the WZA exhibited strong performance

It is always worth keeping in mind that if there is true local adaptation influencing a particular gene, it may generate linkage disequilibrium (LD) with closely linked sites, so multiple genes could exhibit similar GEA signals. Though thus is a general problem with genome scans, not specific to the WZA, top-candidate method or even GEA analyses. GEA analyses are association based analyses, getting to causation from GEA results additional work.

## Note
I was using this repository while developing and performing analyses for the project described in our preprint. As a result a lot of what is contained in this repository is specific to our project exploring the properties of the WZA.
