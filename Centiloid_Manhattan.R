library(tidyverse)
library(qqman)



#Reading in CHR21 and autosome results
ABCDS_Centiloid.assoc <- read.table("~/Documents/Thesis_Data/ABCDS_Centiloid.assoc.frq", header=TRUE, quote="", comment.char="")

MEQTL_Results_Centiloid <- read.delim("~/Documents/Thesis_Data/MEQTL_Results_Centiloid", quote="")

#Using regular expressions to make CHR and BP columns from SNP naming
MEQTL_Results_Centiloid$CHR <- sub(":.*", "", MEQTL_Results_Centiloid$SNP)

MEQTL_Results_Centiloid$BP <- sub(".*:(\\d+)\\[.*", "\\1", MEQTL_Results_Centiloid$SNP)  

#reordering columns to match autosomes, making columns match. 

MEQTL_Results_Centiloid <- MEQTL_Results_Centiloid %>% select(CHR, BP, everything())

MEQTL_Results_Centiloid <- MEQTL_Results_Centiloid %>% rename(P = p.value)

MEQTL_Results_Centiloid <- MEQTL_Results_Centiloid %>% select(CHR,BP,SNP,P)

ABCDS_Centiloid.assoc <- ABCDS_Centiloid.assoc %>% select(CHR,BP,SNP,P)

MEQTL_Results_Centiloid$CHR <- as.numeric(MEQTL_Results_Centiloid$CHR)

MEQTL_Results_Centiloid$BP <- as.numeric(MEQTL_Results_Centiloid$BP)

Centiloid_Complete <- bind_rows(MEQTL_Results_Centiloid, ABCDS_Centiloid.assoc)

Centiloid_Complete <- Centiloid_Complete %>% arrange(CHR,BP)

#QQ plotting
qq(Centiloid_Complete$P, main= "Centiloid") 


#Removing uneccesary values for quicker plotting 
Centiloid_Complete <- Centiloid_Complete %>% 
  filter(-log10(P)>1)



manhattan(Centiloid_Complete, chr="CHR", bp="BP", snp="SNP", p="P", main= "Centiloid", annotatePval = 5e-8)

