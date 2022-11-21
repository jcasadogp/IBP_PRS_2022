# Load different types of libraries.
library(remotes)
library(data.table)
library(magrittr)
library(bigsnpr)
library(fmsb)
library(runonce)

# Read in all the different files.
sumstats <- fread(snakemake@input[[1]])
corr <- fread(snakemake@input[[2]])
pcs <- fread(snakemake@input[[3]])
phenotype <- fread(snakemake@input[[4]])
bed <- snp_readBed(snakemake@input[[5]])
obj.bigSNP <- snp_attach(bed)
info <- readRDS(snakemake@input[[6]])
bfile <- snakemake@params[[1]]


corr$FID <- as.character(corr$FID)
corr$IID <- as.character(corr$IID)

# Filter out hapmaps
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

# Create necessary variables
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
ld <- "EUR.hg38"
fam.order <- NULL

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
# calculate LD
 #Get maximum amount of cores
#NCORES <- nb_cores()
options(bigstatsr.check.parallel.blas = FALSE)
for (chr in 1:22) {
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = nb_cores(),
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- phenotype[fam.order, on = c("FID", "IID")]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Height~Sex+", .) %>%
    as.formula %>%
    glm(., data = y, family=binomial) %>%
    summary
null.r2 <- fmsb::NagelkerkeR2(null.model)


beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
prs_results <- big_prodVec(    genotype,
                            beta_inf,
                            ind.row = ind.test,
                            ind.col = info_snp$`_NUM_ID_`)

reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Height~PRS+Sex+", .) %>%
    as.formula
reg.dat <- y
reg.dat$PRS <- pred_inf
inf.model <- lm(reg.formula, dat=reg.dat) %>%
    summary
(result <- data.table(
    infinitesimal = inf.model$r.squared - null.r2,
    null = null.r2
))

# Creating the output. 
prs_results <- reg.dat$PRS
max_r_sqr <- inf.model
prs_results <- NULL
max_r_sqr <- NULL
write.table(prs_results, quote = FALSE,snakemake@output[[1]])
write.table(max_r_sqr, quote = FALSE,snakemake@output[[2]])




