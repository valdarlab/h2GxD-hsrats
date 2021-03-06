# h2GxD-hsrats

This document describes how to reproduce heritability analysis of various phenotypes in response to high- or low-fat diet in heterogeneous stock rats, as performed by Deal et al (2022). 

Sex and genetic specific effects on behavioral, but not metabolic, responses to a high fat diet in heterogeneous stock rats   
========================================================================

Environment prep
----------------

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
Run the following R script to transform and clean up the raw data located in the `source_data` directory: 

```
Rscript data_transformations.R
```

The following files (one for each sex) will be produced in the `derived_data` directory: `transformed_data_female-G3package.csv` and `transformed_data_male-G3package.csv`. 

Run the following Rscript to create additive relationship matrices:

```
Rscript relationship_matrices.R 
```

Analysis
--------

Run the following bash script to build the models and calculate model statistics. For each sex/phenotype combination, this script will spin off a background job that calls `pheno_heritability_est.R`. 

```
bash heritability_est.sh 
```

The default behavior is to create background jobs using standard bash commands. To create the job via SLURM, set the `-m` flag to "slurm": 

```
bash heritability_est.sh -m slurm 
```

In either case, a `logs` directory will be created where stdout from each job will print to its own `.log` file. You can track the progress of the program for each phenotype here. 

When analysis is complete, flat files will be produced in the `posterior_samples` directory and csv files with summary statistics will be produced in the `heritability_stats` directory for each sex/phenotype combination.

Results
-------

Create a results directory, and combine the csv files from the previous step into one file: 

```
mkdir results
awk 'FNR==1 && NR!=1 {next;}{print}' derived_data/heritability_stats/*.csv > results/combined.csv
```

Generate a report of Bayes-Factor analysis:

```
Rscript -e 'library(rmarkdown); rmarkdown::render("SavageDickeyGxT.Rmd", "pdf_document")'
```

In addition to the `SavageDicketGxT.pdf` report, `results/BF_results.csv` will be produced and `results/plots/` and `results/posterior_densities/` directories will be populated with one file per sex/phenotype combination. 

Notes
-----
These steps must be performed in the order provided to preserve dependencies. 

If at any point you want to open an RStudio session from within the Docker container, go to [http://localhost:8787/](http://localhost:8787/) and login with user `rstudio` and password `pw123`. 