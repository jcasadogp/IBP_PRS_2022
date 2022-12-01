# IBP_PRS_2022

This project aims to compare four different tools that calculate Polygenic Risk Scores (PRS). These tools are PLINK, PRSice, Lasso Sum and LDPred. Currently, this workflow only allows for testing of **binary** phenotypes, such as cases and controls.

The heart of this project is a Snakemake script, therefore allowing for extension with more tools. Additionally, the project is supported for users that have access to a remote cluster. The scripts require the submission of a SLURM job, however if your cluster uses another job scheduler this can easily be adapted as well. The aim of using this pipeline is not to take the PRS directly. Rather it should aid the user in selecting the best tool to use for calculating PRS on their specific dataset. The user should rerun their own analysis on such tool with their data to allow for more fine tuning. 

_If you are working in group directory, remember to give access to the other members after creating any file or directory._

## PREVIOUS REQUIREMENTS

* Conda

## DOWNLOAD THIS REPOSITORY

In order to run this tool you have to download this GitHub repository.
```
git clone https://github.com/jcasadogp/IBP_PRS_2022.git
```

#### ERRORS or WARNINGS I GET WHEN I TRY ALL THESE COMMANDS:
![Captura de pantalla 2022-12-01 a las 12 30 01](https://user-images.githubusercontent.com/80517901/205041802-03e42839-e716-4641-a418-c7e6643ef194.png)

* We need to be careful with case-sensitive paths.
 
![Captura de pantalla 2022-12-01 a las 12 34 17](https://user-images.githubusercontent.com/80517901/205042720-ba0a9b6f-1aa9-4590-9637-c694f0f46711.png)

* The command is NOT CLONE, because I have been able to edit a push changes from local, and we don't want that.

It includes the following folders:
* Conda environments: several .yml files that will be used by the different scripts to activate the contained conda environments. These .yml files are what allows for the running of lassoSum, LDpred and generation of performance metric plots.
* plink: installation and executable files for PLINK
* prsice: installation, executable and R script for PRSice
* r_scripts: it contains several R scripts that are used in the project 
* snakemake: it contains the global Snakefile that will run the four tools as well as the job file (.pbs) that must be run by the user.


## INPUT DATA:

Data needed:

* GWAS summary statistics
* 1000 genome file to use as reference, which should be a subset of the cohort that is of the same ancestral background as your target cohort
* Target files: PLINK binary format files (.fam, .bed, .bim)
* Phenotype file: File containing the phenotype file with FID, IID, and as many columns as phenotypes you want to look at.
* Covariates and eigenvectos files (eigenvalues are not really necessary):

All this data will be placed in the data/ directory. To create it run:

```
cd IBP_PRS_2022
mkdir data/
cp FILE data/
```

## OUTPUT DATA:

The output data will be placed in a directory called ```output_data/``` that contains a different directry per tool. Also, inside each of the tools' directories, we will separate the data into ```target_data/``` and ```exteranl_data/```, depending if it comes from the tool with the 1000genome data (```external_data/```) or from the target data (```target_data/```). To create all these directories run:

```
cd IBP_PRS_2022
mkdir output_data/
cd output_data/
mkdir 001_plink/
mkdir 001_plink/target_data/
mkdir 001_plink/external_data/
mkdir 002_prsice/
mkdir 002_prsice/target_data/
mkdir 002_prsice/external_data/
mkdir 003_lassoSum/
mkdir 003_lassoSum/target_data/
mkdir 003_lassoSum/external_data/
mkdir 004_LDpred/
mkdir 004_LDpred/target_data/
mkdir 004_LDpred/external_data/
mkdir 005_comparison/
```


## Adapting the Snakemake script
Snakemake is a wonderful workflow engine which allows for easy adaptability and extension. At the top of the script you will see several global variables, all of which you will need to fill in. They include:

```
# === Prefix of the files  ===
TARGET_files = "IBD_GSA_fin"
EXTERNAL_file = "1000G_EUR_fin"
GWAS_file = "GWAS_summary_stats.txt"
PHENO_file = "final_phenotypes.txt"

# === Parameters for tools ===
lasso_thresholding_values = "0.2,0.5,0.9,1" #must be in a comma separated format with values from 0 to 1
genome_build = "EUR.hg38" #options: "EUR.hg19", "AFR.hg19", "ASN.hg19", "EUR.hg38", "AFR.hg38", "ASN.hg38"
```

The user should provide the file names respectively. The parameter arguments can be specified for user's specific target datas genome build. Further, the thresholding and shrinkage parameters can be adapted to allow for more niche testing.

## Adapting the job script




# GWAS FILE -> should we tell the user to recode the column names?


## Output


### Files


### Plots
Several plots will be generated for the prediction of PRS and comparison of the performance metrics across tools.
* Each tool will have a boxplot of PRS generated comparing the cases and controls
* A ROC curve will be generated to show how well the logistic model built from each tools PRS predictions is at correctly classifying cases and controls
* Bar plots comparing the AUC and R<sup>2</sup> values across the tools, generated using 1000 bootstrap samples

### Sources

Berisa, T. & Pickrell, J. K. Approximately independent linkage disequilibrium blocks in human populations. Bioinformatics 32, 283-285 (2015).

Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020). 

Mak, Timothy Shin Heng, et al. "Polygenic scores via penalized regression on summary statistics." Genetic epidemiology 41.6 (2017): 469-480. 
* [lassoSum](https://github.com/tshmak/lassosum) repository


