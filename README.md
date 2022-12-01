# IBP_PRS_2022

This project aims to compare four different tools that calculate Polygenic Risk Scores (PRS).

If you are working in group, remember to give access to the other members after creating any file or directory.

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
* Conda environments: several .yml files that will beused by the different scripts to create conda environments. (JULIA: The user needs conda right??)
* plink: installation and executable files for PLINK
* prsice: installation, executable and R script for PRSice
* r_scripts: it contains several R scripts that are used in the project 
* snakemake: it contains the global Snakefile that will run the four tools as well as the job file (.pbs) that must be run by the user.


## INPUT DATA:

Data needed:

* GWAS summary statistics
* 1000 genome file to use as reference
* Target files: PLINK format files (.fam, .bed, .bim)
* Phenotype file: File containing the phenotype file with FID, IID, and as many columns as phenotypes you want to look at.
* Covariates and eigenvectos files (eigenvalues are not really necessary):

All this data will be placed in the data/ directory. To create it run:

```
cd DIRECTORY DOWNLOADED WITH THE PARENT DIRECTORY
mkdir data/
cp FILE data/
```

## OUTPUT DATA:

The output data will be placed in a directory called output_data/ that contains a different directry per tool. Also, inside each of the tools' directories, we will separate the data into target_data/ and exteranl_data/, depending if it comes from the tool with the 1000genome data (external_data/) or from the target data (target_data/). To create all these directories run:

```
cd DIRECTORY DOWNLOADED WITH THE PARENT DIRECTORY
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




