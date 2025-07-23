### Mutation rates in coding vs. non-coding regions
### AVC Winter 2024

#### Set up environment --------------------------------------------------------

library(tidyverse)

# get metadata
genome_metadata <- read.csv("genome_metadata.csv") %>% 
  filter(reference == "Kbre")

#### Get Pmac chromosome names -------------------------------------------------

kbre_chrom <- read.csv("SNPs_by_chromosome_Kbre.csv") %>% 
  dplyr::select(GenBank, chrom.num) %>% 
  distinct()

#### Create BED file for coding regions ----------------------------------------

pmac_code <- read.delim("G:/My Drive/00 UW/00.5 W.A.D.E. lab resources/Intern Projects/Mutation rate variability/06 coding_noncoding_snps/Pmac_protein_coding_regions.tsv", sep = "\t") %>% 
  dplyr::select(Chromosome, Begin, End, Name) %>% 
  filter(Chromosome != "MT") %>% 
  filter(Chromosome != "") %>% 
  left_join(pmac_chrom, by = c("Chromosome" = "Chromosome.name")) %>% 
  dplyr::select(-Chromosome) %>% 
  rename("Chromosome" = RefSeq.seq.accession) %>% 
  relocate(Chromosome, .before = Begin)

kbre_genes <- read.delim("Kbre_annotations.tsv", sep = "\t") %>% 
  dplyr::select(Chromosome, Begin, End, Name, Gene.Type) %>% 
  filter(Chromosome != "MT") %>% 
  filter(Chromosome != "") %>% 
  left_join(kbre_chrom, by = c("Chromosome" = "chrom.num")) %>% 
  dplyr::select(-Chromosome) %>% 
  rename("Chromosome" = GenBank) %>% 
  relocate(Chromosome, .before = Begin) 

kbre_code <- kbre_genes %>% 
  filter(Gene.Type == "protein-coding")

write_delim(kbre_code, file = "Kbre_protein_coding_regions.bed", delim = "\t", col_names = FALSE,
            quote = "none")

# calculate total length of coding regions
code_length <- kbre_code %>% 
  mutate(length = End - Begin) %>% 
  summarize(tot.len = sum(length)) %>% 
  pull(tot.len)

#### Create BED file for non-coding regions ------------------------------------

kbre_noncode <- kbre_genes %>% 
  filter(Gene.Type != "protein-coding")

write_delim(kbre_noncode, file = "Kbre_non_coding_regions.bed", delim = "\t", col_names = FALSE,
            quote = "none")

# calculate total length of non-coding regions
noncode_length <- kbre_noncode %>% 
  mutate(length = End - Begin) %>% 
  summarize(tot.len = sum(length)) %>% 
  pull(tot.len)

################################################################################
#### From here, use the bed files to sort SNPs into coding and non-coding   ####
#### regions using bedtools, and count the number of SNPs in both regions   ####
#### using awk. Bash code for both of these steps can be found in           ####
#### code_noncode_readdepth_miscbash.qmd. Once these steps are completed,   ####
#### the code below can be used to read in the txt files, calculate         ####
#### mutation rate, and plot mutation rate in coding and non-coding regions.####
################################################################################

#### Read in txt files with number of SNPs in coding and noncoding regions -----

coding_data_files <- list.files(path = "M:/Mutation_Rate/snp_counts",
                                pattern = "_coding_snps.txt")

pcoding_data <- do.call(rbind, lapply(paste0("M:/Mutation_Rate/snp_counts/", coding_data_files), 
                                     read.delim, sep=" ", 
                                     header = FALSE, 
                                     col.names = c("hets","homs","aallele"))) %>% 
  filter(hets > 0) %>% 
  mutate(file = coding_data_files) %>% 
  separate(file, into = c("sp_code", NA), sep = "_") %>% 
  mutate(tot.length = code_length) %>% 
  mutate(snp_type = "coding")

noncoding_data_files <- list.files(path = "M:/Mutation_Rate/snp_counts",
                                   pattern = "noncoding_snps.txt")

noncoding_data <- do.call(rbind, lapply(paste0("M:/Mutation_Rate/snp_counts/", noncoding_data_files), 
                                      read.delim, sep=" ", 
                                      header = FALSE, 
                                      col.names = c("hets","homs","aallele"))) %>% 
  filter(hets >0) %>% 
  mutate(file = noncoding_data_files) %>% 
  separate(file, into = c("sp_code", NA), sep = "_") %>% 
  mutate(tot.length = noncode_length) %>% 
  mutate(snp_type = "noncoding")

combined_data <- bind_rows(pcoding_data, noncoding_data) %>% 
  left_join(genome_metadata, by = c("sp_code" = "abbrev"))

#### Calculate mutation rate for coding and non-coding regions -----------------

code_noncode_data <- data.frame()

for (i in 1:nrow(combined_data)) {
  AALLELE=combined_data$aallele[i]
  TOTAL=combined_data$tot.length[i] - AALLELE
  HET=combined_data$hets[i]
  HOMALT=combined_data$homs[i]
  HOMREF=TOTAL-(HOMALT+HET)
  
  # Heterozygosity
  PI=HET/TOTAL
  
  # Divergence (number of alt alleles in hets and hom. alt divided by the total number of 
  # alleles)
  DIV=(HET+(2*HOMALT))/(2*TOTAL)
  
  #Time to coalesence in years
  T=combined_data$div_mcgowen2020[i]
  
  #noPi mutation rate
  noPI_mutation <- DIV/(2*T)
  
  #Pi mutation rate
  PI_mutation <- (DIV-PI)/(2*T)
  
  temp_coding_data <- data.frame(species = rep(combined_data$sp_code[i],2), 
                                snp_type = rep(combined_data$snp_type[i],2), 
                                method = c("noPI","PI"), 
                                rate = c(noPI_mutation, PI_mutation))
  
  code_noncode_data <- rbind(code_noncode_data, temp_coding_data)
}

code_noncode_data <- code_noncode_data %>% 
  mutate(suborder = case_when(species %in% c("Bmus","Egla","Bacu","Erob","Bric", "Mnov")~"Mysticete",
                   species %in% c("Kbre","Igeo","Pele","Mden","Oorc","Psin","Ddel",
                                  "Hamp","Scoe","Ggri","Gmel","Lalb")~"Odontocete",
                   TRUE~NA)) %>% 
  filter(species != "Kbre")

#### Test for differentiation in rates in coding and non-coding regions --------

coding_normality <- shapiro.test(code_noncode_data$rate)
coding_normality

code_noncode_aov <- aov(rate ~ snp_type, data = code_noncode_data %>% 
                          filter(method == "noPI"))
summary(code_noncode_aov)

eta_squared(code_noncode_aov, partial = FALSE)

order_code_noncode_aov <- code_noncode_data %>% 
  filter(method == "noPI") %>% 
  group_by(suborder) %>% 
  tidyr::nest() %>%
  dplyr::mutate(.data = .,
                aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(rate ~ snp_type, data = .))))

myst_coding_aov <- order_code_noncode_aov$aov_results[[1]]
odon_coding_aov <- order_code_noncode_aov$aov_results[[2]]

coding_difference <- code_noncode_data %>% 
  filter(method == "noPI") %>% 
  pivot_wider(names_from = "snp_type", values_from = "rate") %>% 
  mutate(diff_rate = coding-noncoding)
  
#### graph time! ---------------------------------------------------------------
library(ggplot2)
library(PNWColors)

# data frame for previous estimations of mutation rate
est_mr <- data.frame(yint=1.08e-8)

ggplot(data=code_noncode_data, aes(x=snp_type,y=rate,shape=method,color=species)) +
  geom_point(size=3.5) +
  theme_light() +
  labs(x="Coding vs. noncoding regions", y="Mutations/site/generation", 
       title="Mutation Rate by region") +
  scale_shape_discrete(labels=c('w/o ancestral\nheterozygosity', 'w/ ancestral\nheterozygosity'), name="Method") +
  scale_color_manual(values=pnw_palette(n=18,name="Sunset2"), 
                     #labels=c("Risso's dolphin","Short-finned\npilot whale","Indo-Pacific\nhumpback dolphin"),
                     name="Species") +
  theme(axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=5,l=0))) +
  theme(text = element_text(size = 20)) +
  geom_hline(data=est_mr, mapping=aes(yintercept=yint, linetype="A")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

# optional violin plot, data aggregated across all species
ggplot(data=code_noncode_data, aes(x=suborder,y=rate,fill=snp_type)) + 
  geom_violin() +
  theme_light() +
  labs(x="Coding vs. noncoding regions", y="Mutations/site/generation", 
       title="Mutation Rate in coding vs. noncoding regions") +
  #scale_shape_discrete(labels=c('w/o ancestral\nheterozygosity', 'w/ ancestral\nheterozygosity'), name="Method") +
  scale_fill_manual(values=pnw_palette(n=3,name="Sunset2"), 
                     #labels=c("Risso's dolphin","Short-finned\npilot whale","Indo-Pacific\nhumpback dolphin"),
                     name="SNP locus") +
  theme(axis.title.y = element_text(margin = margin(t=0,r=8,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=8,r=0,b=5,l=0))) +
  theme(text = element_text(size = 20))

save(coding_difference, code_noncode_aov, 
     code_noncode_data, odon_coding_aov, 
     myst_coding_aov, file = "Mutation-rate-coding-noncoding.Rdata")
