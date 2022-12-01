# IBP_PRS_2022

This project aims to compare four different tools that calculate Polygenic Risk Scores (PRS).

If you are working in group, remember to give access to the other members after creating any file or directory.

## INPUT DATA:

Data needed:

* GWAS summary statistics
* 1000 genome file to use as reference
* Target files: PLINK format files (.fam, .bed, .bim)
* Phenotype file: File containing the phenotype file with FID, IID, and as many columns as phenotypes you want to look at.
* Covariates and eigenvectos files (eigenvalues are not really necessary):

All this data will be placed in the data/ directory. To create it run:

```
mkdir data/
cp FILE data/
```

## OUTPUT DATA:

The output data will be placed in a directory called output_data/ that contains a different directry per tool. Also, inside each of the tools' directories, we will separate the data into target_data/ and exteranl_data/, depending if it comes from the tool with the 1000genome data (external_data/) or from the target data (target_data/). To create all these directories run:

```
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




