# Chromosome level mutation rate estimation with cetaceans
# Sophie Garrote and Amy Van Cise
# 28 Feb 2024

# mutation rate calculation script derived from Robinson et al. 2022

# Counting/processing of chromosome-level SNPs done in bash

################################################################################

library(tidyverse)

setwd("M:/Mutation_Rate")

### Get Pmac chromosome length -------------------------------------------------

#create length column
chromfile <- read.delim("Vcf_mutation/Pmac_chrom_total.txt", col.names = "filename", header = FALSE)

chrom.length <- chromfile %>%
  separate(filename, into = c(NA, "length"), sep = ",") %>%
  separate(length, into = c(NA, "length"), sep = "=") %>%
  mutate(length = parse_number(length)) %>% 
  mutate(chrom.num = row_number()) %>% 
  mutate(chrom.num = case_when(chrom.num == max(chrom.num) ~ "X", TRUE ~ as.character(chrom.num)))

### Get genome metadata --------------------------------------------------------

genome_metadata <- read.csv("R_script_info_mrate.csv")

### Create data frame from all text files --------------------------------------

# accessing processed files containing SNPs partitioned by chromosome
data_files <- list.files("Vcf_mutation", pattern = "chrom_count.txt")

# Read the files in, add Pmac chrom.length
data_files_list <- lapply(paste0("Vcf_mutation/",data_files), 
                          function(x) {read.delim(file = x, header = FALSE, sep =" ") %>% 
                          dplyr::rename("abbrev" = 1, "reference" = 2, "hets" = 3, 
                          "hom" = 4, "aallele" = 5)} %>% cbind(chrom.length))

# Combine them
chrom_snp_data <- do.call("rbind", lapply(data_files_list, as.data.frame)) 

# add info for mutation rate calculations from csv file
chrom_snp_data <- chrom_snp_data %>% 
  left_join(genome_metadata, by = c("abbrev", "reference"))

#data frame for chromosome mutation rates---------------------------------------
species_chrom_mutation_data <- data.frame()

### loop that iterates through mutation rate calculation w/ info from each row 
  for (i in 1:nrow(chrom_snp_data)) {
    AALLELE=chrom_snp_data$aallele[i]
    TOTAL=chrom_snp_data$length[i] - AALLELE
    HET=chrom_snp_data$hets[i]
    HOMALT=chrom_snp_data$hom[i]
    HOMREF=TOTAL-(HOMALT+HET)
    
    # Heterozygosity
    PI=HET/TOTAL
    
    # Divergence (number of alt alleles in hets and hom. alt divided by the total number of 
    # alleles)
    DIV=(HET+(2*HOMALT))/(2*TOTAL)
    
    #Time to coalesence
    T=chrom_snp_data$div_timetree[i]/chrom_snp_data$gen_time[i]
    
    #noPi mutation rate
    noPI_mutation <- DIV/(2*T)
    
    #Pi mutation rate
    PI_mutation <- (DIV-PI)/(2*T)
    
    temp_chrom_data <- data.frame(species = rep(chrom_snp_data$abbrev[i],2), 
                                  chrom = rep(chrom_snp_data$chrom.num[i],2), 
                                  method = c("noPI","PI"), 
                                  rate = c(noPI_mutation, PI_mutation))
                                  
    species_chrom_mutation_data <- rbind(species_chrom_mutation_data, temp_chrom_data)
  }

################################################################################

#changing character vector to factor to correct graph chrom order
species_chrom_mutation_data$chrom <- factor(species_chrom_mutation_data$chrom, 
                                            levels=c("1","2","3","4","5","6","7","8","9","10","11",
                                                     "12","13","14","15","16","17","18","19","20","X"))
## graph time!
library(ggplot2)
library(PNWColors)
install.packages("viridis")
library(viridis)

# data frame for previous estimations of mutation rate
est_mr <- data.frame(yint=1.08e-8)

ggplot(data=species_chrom_mutation_data, aes(x=chrom,y=rate,shape=method,color=species)) +
  geom_point(size=3.5) +
  theme_light() +
  labs(x="Chromosome", y="Mutations/site/generation", 
       title="Chromosome-Level Mutation Rate") +
  scale_shape_discrete(labels=c('w/o ancestral\nheterozygosity', 'w/ ancestral\nheterozygosity'), name="Method") +
  scale_color_viridis(discrete = TRUE) +
  theme(axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=5,l=0))) +
  theme(text = element_text(size = 20)) +
  geom_hline(data=est_mr, mapping=aes(yintercept=yint, linetype="A")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

# OPTIONAL:what is the range in mutation rate for species vs. chromosomes?
  species_chrom_mutation_data %>% 
    group_by(species) %>% 
    summarise(minmax = range(rate)) %>% 
    group_by(species) %>% 
    summarise(range = max(minmax)-min(minmax))
  
  species_chrom_mutation_data %>% 
    group_by(chrom) %>% 
    summarise(minmax = range(rate)) %>% 
    group_by(chrom) %>% 
    summarise(range = max(minmax)-min(minmax))
  
  