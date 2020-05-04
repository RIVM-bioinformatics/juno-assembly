## Assembly pipeline

W.I.P

To derive valuable information and quality metrics from the paired end Illumina data, the data are analyzed in a sequential order using various algorithms to assemble the genome and to decide which samples can participate in analysis of antimicrobial resistance genes, virulence genes and perform cluster analyses and which samples are excluded from further analysis due low quality or contaminated isolates.

First FastQC (Andrews, et al.,2010) is performed to assess the quality of the raw Illumina reads. Trimmomatic (Bolger, et al., 2014) is used to remove poor quality data and adapter sequences. The Trimmomatic parameters used are a sliding window of 5:30 and a minlen config value of 50. FastQC is used once more to assess the quality of the trimmed reads. And Picard (https://broadinstitute.github.io/picard/) determines the library fragment lengths.

These reads are assembled into scaffolds by SPAdes (Nurk, et al., 2017) by means of de novo assembly of the genome. For de novo assembly the SPAdes Isolate option was used and Kmer sizes of 21,33,55,77 and 99 were used. The scaffolds that are lower than the threshold of 500 nucleotides are filtered out. Scaffolds with the length of 500 or more continue in the upstream analysis. 

QUAST (Gurevich, et al., 2013) is performed to assess the quality of the filtered scaffolds. To assess the quality of the microbial genomes recovered from isolate CheckM (Parks, et al., 2014) is used. CheckM gives scores for completeness, contamination and strain heterogeneity. Bbtools (Bushnell, 2014) is performed to generate scaffold alignment metrics. 

As the last step of the pipeline MultiQC (Ewels, et al., 2016) is used to fit and summarize analysis results and quality assessments in a single report for dynamic visualization.



```

Dry run:
```
bash start_here.sh -n
```

Start pipeline:
```
bash start_here.sh
```
