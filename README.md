# h2GxD-hsrats

This document describes how to reproduce heritability analysis of various phenotypes in response to high- or low-fat diet in heterogeneous stock rats, as performed by Deal et al (2022). 

Sex and genetic specific effects on behavioral, but not metabolic, responses to a high fat diet in heterogeneous stock rats   
========================================================================

Environment prep
----------------

This project uses a Docker container to produce an environment similar to that used in the original analysis (e.g. R v4.2.1 and associated R package versions). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

Build the docker container:

```
docker build . -t gxd 
```

Run the docker container, opening a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it gxd /bin/bash
```

Navigate to the working directory: 

```
cd home/rstudio/work 
```

Data prep
---------
Run the following R script to summarize and transform/clean up the raw data located in the `source_data` directory: 

```
Rscript data_transformations.R
```

This script will calculate descriptive statistics )mean and standard deviation) for raw phenotypes and save them in `results/phenotype_summary.csv`.

This script will also transform the raw data as specified in `source_data/phenotype_transformations.csv`. The following files (one for each sex, and one combined) will be produced in the `derived_data` directory: `transformed_data_female-G3package.csv`, `transformed_data_male-G3package.csv`, and `transformed_data_all-G3package.csv`. 

Run the following Rscript to create additive relationship matrices:

```
Rscript kinship_using_QTLRel_gen35_microbiome.R 
```

Three pedigrees and three relatedness matrices (male, female and combined) will be created and stored in .txt files in the `derived_data` directory. 

Analysis
--------

Run the following script to perform factor analysis on transformed phenotypes: 

```
Rscript factor_analysis.R
```

This script will perform factor analysis on the transformed dataset with missing values imputed, save the factor loadings to the `results` folder, and calculate factor phenotypes for each factor (the sum of the phenotypic values multiplied by their corresponding factor loadings). The transformed data files in the `derived_data` directory will be updated with these factor phenotypes, and downstream analyses will be performed on these factor phenotypes as well as the original phenotypes. 

Run the following Rscript to perform covariate analysis:

```
Rscript anova.R
```

-------
This script performs two primary analysis:
1. Sex-stratified analysis asks whether there is a significant effect of diet, stratified by sex (**RESULTS FILENAME**) 
2. Covariate analysis. Results of this analysis will be saved in the `results` directory. This includes p-values for the significance of each term in predicting each phenotype, as well as heritability estimates for each phenotype, calculated via ICC. 
---------

Run the following bash script to build the models and calculate model statistics. For each sex/phenotype combination, this script will spin off a background job that calls `pheno_heritability_est.R`. 

```
bash heritability_est.sh 
```

The default behavior is to create background jobs using standard bash commands. To create the job on a cluster via slurm, set the `-m` flag to "slurm": 

```
bash heritability_est.sh -m slurm 
```

The time limit for slurm jobs is specified at 3 days. The longest running phenotype on UNC's high-performance computing cluster was 48 hours. 

In either case, a `logs` directory will be created where stdout from each job will print to its own `.log` file. You can track the progress of the program for each phenotype here. 

When analysis is complete, flat files will be produced in the `posterior_samples` directory and csv files with summary statistics will be produced in the `heritability_stats` directory for each sex/phenotype combination.

Results
-------

Combine the csv files from the previous step into one file, which will be stored in the results directory created by the ANOVA analysis:

```
awk 'FNR==1 && NR!=1 {next;}{print}' derived_data/heritability_stats/*.csv > results/combined.csv
```

Generate a report of Bayes-Factor analysis:

```
Rscript -e 'library(rmarkdown); rmarkdown::render("SavageDickeyGxT.Rmd", "pdf_document")'
```

In addition to the `SavageDicketGxT.pdf` report, `results/BF_results.csv` will be produced and `results/plots/` and `results/posterior_densities/` directories will be populated with one file per sex/phenotype combination. 

**Note** that this step requires almost 15 GB of memory. If you are using Docker for Mac, you will need to increase the memory limit above the default 2 GB.

Notes
-----
These steps must be performed in the order provided to preserve dependencies. 

If at any point you want to open an RStudio session from within the Docker container, go to [http://localhost:8787/](http://localhost:8787/) and login with user `rstudio` and password `pw123`. 