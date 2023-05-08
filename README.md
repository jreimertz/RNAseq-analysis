# RNAseq-analysis

A web application developed in RShiny for performing a number of different 
bioinformatics analyses on an inputted file containing RNAseq data. The 
requirements to use this application include: a tab or comma delimited sample 
information meta data file, a comma or tab delimited normalized counts data
file. Pre-computed differential expression analysis results are optional but
can also be uploaded to the app. <br><br>


## How to use this application

This application can be used to explore the sample data including: investigating
the contents and data types of each column, and generating plots using the 
discrete and continuous data stored in this file. Additionally, users can also
explore the normalized bulk RNA-seq counts matrix including:
applying filters to the data based on the percentile of variance between genes 
and the number of samples with non-zero counts, creating diagnostic plots of
the full counts matrix with the distribution of genes passing vs not passing
filters displayed, and perform PCA with the top N PC's displayed as a 
beeswarm plot. Differential analysis can be performed on the filtered counts
matrix or a results file from previously computed differential analysis can
be uploaded to explore the results and generate a volcano plot displaying the
distribution of log2 fold change and -log10(adjusted p-value). Finally, pairwise
pearson correlation network analysis can be performed on the normalized counts 
data including displaying the visualization of the network and associated 
metrics. <br><br>

## Access to the app
This application can be run locally by downloading the project on this github
repository and has also been published at 
<a href="https://jreimertz.shinyapps.io/main/">RNA-seq Analysis</a>
Additionally, a demo recording of this application use can be found
<a href="https://youtu.be/IfnhUCnUapk"></a><br><br>

## Purpose

This app was originally built as part of the BF591 final project at 
Boston University using the data accessed via 
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810">GSE64810</a>
<br><br>

