library(tidyverse)
library(MatrixEQTL)

#Filtering phenotype data to Centiloid
ABCDS_Centiloid <- ABCDS_FILTERED %>% select(FID,IID,Centiloid)

ABCDS_Centiloid$Sample <- paste(ABCDS_Centiloid$FID, ABCDS_Centiloid$IID, sep = "_")

#Making the Samples appear the same as in the SNP data set for easy comparison
ABCDS_Centiloid <- ABCDS_Centiloid %>% select(Sample, Centiloid)

ABCDS_Centiloid$Sample <- paste0("X", ABCDS_Centiloid$Sample)

#Making centiloid sample wide 

ABCDS_Centiloid <- pivot_wider(ABCDS_Centiloid, names_from = Sample, values_from = Centiloid) 



#Finding common columns 
Centiloid_Col <- colnames(ABCDS_Centiloid)
Snp_Col <- colnames(ABCDS_Snp)

common_columns <- intersect(Cent_Col, Snp_Col)
#Filtering by common columns

ABCDS_Centiloid <- ABCDS_Centiloid %>% select(all_of(common_columns))

ABCDS_Centiloid_Snp <- ABCDS_Snp %>% select(all_of(common_columns))

rownames(ABCDS_Centiloid_Snp) <- ABCDS_Snp$snpID

#Doing the same with covariates
ABCDS_Centiloid_Cov <- ABCDS_FILTERED %>% select(FID,IID,Sex, age_at_visit, PC1,PC2,PC3,PC4)

ABCDS_Centiloid_Cov$Sample <- paste(ABCDS_Centiloid_Cov$FID, ABCDS_Centiloid_Cov$IID, sep = "_")

ABCDS_Centiloid_Cov <- ABCDS_Centiloid_Cov %>% select(Sample,Sex,age_at_visit,PC1,PC2,PC3,PC4)

ABCDS_Centiloid_Cov$Sample <- paste0("X", ABCDS_Centiloid_Cov$Sample)

ABCDS_Centiloid_Cov <- pivot_longer(ABCDS_Centiloid_Cov, cols = -Sample) %>% pivot_wider(names_from = Sample)

ABCDS_Centiloid_Cov <- ABCDS_Centiloid_Cov %>% select(all_of(common_columns))


rownames(ABCDS_Centiloid) <- "Centiloid"

ABCDS_Centiloid <- as.matrix(ABCDS_Centiloid)
ABCDS_Centiloid_Cov <- as.matrix(ABCDS_Centiloid_Cov)
ABCDS_Centiloid_Snp <- as.matrix(ABCDS_Centiloid_Snp)



write.table(ABCDS_Centiloid, "ABCDS_Centiloid_Geno.txt", sep = "\t", quote = F)
write.table(ABCDS_Centiloid_Cov, "ABCDS_Centiloid_Cov.txt", sep = "\t", quote = F)
write.table(ABCDS_Centiloid_Snp, "ABCDS_Centiloid_Snp.txt", sep = "\t", quote = F)


snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA";
snps$LoadFile("ABCDS_Centiloid_Snp.txt")

gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA";
gene$LoadFile("ABCDS_Centiloid_Geno.txt")

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA";
cvrt$LoadFile("ABCDS_Centiloid_Cov.txt")

useModel=modelLINEAR

pvOutputThreshold= 1

errorCovariance = numeric()

Matrix_eQTL_engine(snps = snps, gene = gene, 
                   cvrt = cvrt, output_file_name = "MEQTL_Results_Centiloid", useModel=useModel, 
                   pvOutputThreshold=pvOutputThreshold, errorCovariance=errorCovariance)

