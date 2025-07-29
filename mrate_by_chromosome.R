# Chromosome level mutation rate estimation with cetaceans
# Sophie Garrote and Amy Van Cise
# 28 Feb 2024

# mutation rate calculation script derived from Robinson et al. 2022

# Counting/processing of chromosome-level SNPs done in bash

################################################################################

library(tidyverse)
library(PNWColors)
library(viridis)
library(dunn.test)

### Get Pmac chromosome length -------------------------------------------------

#create length column
chromfile <- read.delim("Kbre_chrom_total.txt", header = TRUE)

chrom.length <- chromfile %>%
  filter(Chromosome != "AJ554055.1") %>% 
  mutate(chrom.num = row_number()) %>% 
  mutate(chrom.num = case_when(chrom.num == max(chrom.num) ~ "X", TRUE ~ as.character(chrom.num))) %>% 
  rename("length" = RefSeq) %>% 
  dplyr::select(Chromosome, GenBank, chrom.num, length)

### Get genome metadata --------------------------------------------------------

genome_metadata <- read.csv("genome_metadata.csv")

### Get chromosome-level SNP data ----------------------------------------------

chrom_snp_data <- read.csv("SNPs_by_chromosome_Kbre.csv") %>% 
  filter(abbrev != "Kbre") %>% 
  mutate(length = parse_number(length))

# # accessing processed files containing SNPs partitioned by chromosome
# data_files <- list.files("M:/Mutation_Rate/snp_counts", pattern = "chrom_snpcount.txt")
# 
# # Read the files in, add Pmac chrom.length
# data_files_list <- lapply(paste0("M:/Mutation_Rate/snp_counts/",data_files),
#                           function(x) {read.delim(file = x, header = FALSE, sep =" ") %>%
#                           dplyr::rename("abbrev" = 1, "reference" = 2, "blank" = 3,
#                           "hets" = 4, "hom" = 5, "aallele" = 6)} %>% 
#                             filter(hom > 0) %>% 
#                             cbind(chrom.length) %>% 
#                             dplyr::select(-blank))
# 
# # Combine them
# chrom_snp_data <- do.call("rbind", lapply(data_files_list, as.data.frame))
# 
# # add info for mutation rate calculations from csv file
# chrom_snp_data <- chrom_snp_data %>%
#   left_join(genome_metadata, by = c("abbrev", "reference"))
# 
# # save as a csv for later use
# write.csv(chrom_snp_data, file = "SNPs_by_chromosome_Kbre.csv", row.names=FALSE)

# Calculate chromosome mutation rates ------------------------------------------
species_chrom_mutation_data <- data.frame()

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
    T=chrom_snp_data$div_mcgowen2020[i]/chrom_snp_data$gen_time[i]
    
    #noPi mutation rate in generations
    noPI_mutation <- DIV/(2*T)
    
    #Pi mutation rate in generations
    PI_mutation <- (DIV-PI)/(2*T)
    
    #noPi mutation rate in years
    noPI_mutation_years <- DIV/(2*chrom_snp_data$div_mcgowen2020[i])
    
    #Pi mutation rate in years
    PI_mutation_years <- (DIV-PI)/(2*chrom_snp_data$div_mcgowen2020[i])
    
    temp_chrom_data <- data.frame(species = chrom_snp_data$abbrev[i], 
                                  chrom = chrom_snp_data$chrom.num[i], 
                                  method = c("noPI","PI","noPI","PI"), 
                                  time_step = c("gen","gen","year","year"),
                                  rate = c(noPI_mutation, PI_mutation,
                                           noPI_mutation_years, PI_mutation_years))
                                  
    species_chrom_mutation_data <- rbind(species_chrom_mutation_data, temp_chrom_data)
  }

################################################################################

#changing character vector to factor to correct graph chrom order
species_chrom_mutation_data$chrom <- factor(species_chrom_mutation_data$chrom, 
                                            levels=c("1","2","3","4","5","6","7","8","9","10","11",
                                                     "12","13","14","15","16","17","18","19","20","X"))
# Filter to include noPi and annual mutation rate only

mrate_by_chrom <- species_chrom_mutation_data %>% 
  filter(method == "noPI") %>% 
  filter(time_step == "year")

# test for normality
chrom_normality <- shapiro.test(mrate_by_chrom$rate)
chrom_normality
#p-value = 2.82e-06
# not normally distributed
hist(mrate_by_chrom$rate)

# KW test for differentiation
kw_chrom <- kruskal.test(rate~chrom, data = mrate_by_chrom)

mrate_by_chrom %>% kruskal_effsize(rate~chrom) %>% pull(effsize)

# post hoc dunn test
chrom_dunn <- data.frame(dunn.test(mrate_by_chrom$rate, mrate_by_chrom$chrom,
          method = "bonferroni")) %>% 
  filter(P.adjusted < 0.05) %>% 
  separate(comparisons, into = c("chrom1", "chrom2"), sep = " - ")
  

chrom_dunn_diff1 <- chrom_dunn %>% 
  group_by(chrom1) %>% 
  summarise(n_sig1 = n())

chrom_dunn_diff <- chrom_dunn %>% 
  group_by(chrom2) %>% 
  summarise(n_sig2 = n()) %>% 
  full_join(chrom_dunn_diff1, by = c("chrom2" = "chrom1")) %>% 
  mutate(across(c(n_sig1,n_sig2), 
                ~ case_when(is.na(.) ~ 0, TRUE ~ .))) %>% 
  mutate(tot_sig = n_sig1 + n_sig2) %>% 
  filter(tot_sig > 10) %>% 
  mutate(y = 3.5e-10)

## graph time!

ggplot(data=mrate_by_chrom) +
  geom_violin(aes(x=chrom,y=rate,color=chrom, fill=chrom), alpha = 0.5) +
  geom_point(data=chrom_dunn_diff, aes(x=chrom2, y = y), shape = 8, size = 1.5, stroke = 1.1, color = "indianred3") +
  theme_light() +
  labs(x="Chromosome", y="Mutations/site/year") +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=5,l=0))) +
  theme(text = element_text(size = 16), legend.position = "none") +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

# Mean range in mutation rate variabilty among chromosomes ---------------------
chrom_rate_range <- mrate_by_chrom %>% 
    group_by(species) %>% 
    mutate(range = max(rate)-min(rate)) %>% 
    ungroup() %>% 
    summarise(mean = mean(range), max = max(range), min = min(range))

# Variability by chromosome size -----------------------------------------------

chromfile <- chromfile %>% 
  rownames_to_column("chrom_num")

chrom_size_rate <- mrate_by_chrom %>% 
  left_join(chromfile, by = c("chrom" = "chrom_num")) %>% 
  rename("chrom_size" = RefSeq) %>% 
  mutate(chrom_size = parse_number(chrom_size)) %>% 
  filter(chrom != "X")

sizerate_model <- lm(rate ~ chrom_size, data = chrom_size_rate)
summary(sizerate_model)

ggplot(data = chrom_size_rate, aes(x = chrom_size, y = rate)) +
  geom_point() +
  geom_smooth(method = "lm")

save(chrom_rate_range, mrate_by_chrom, 
     chrom_normality, chrom_dunn_diff,
     chrom_size_rate, sizerate_model,
     kw_chrom, file = "mrate_by_chromosome.Rdata")

