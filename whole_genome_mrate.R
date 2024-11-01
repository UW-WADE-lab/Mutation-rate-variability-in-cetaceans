# Calculate whole-genome mutation rate
# AVC 07/2024, edited from Sophe Garrottee 11/2023

###############################################################################

library(tidyverse)
library(ggpubr)
library(PNWColors)
library(png)
library(dunn.test)

#### Load whale images ---------------------------------------------------------
fin <- readPNG("fin-whale-silhouette.png")

#### Load genome metadata, reference metadata, and cetacean heterozygosity -----

genome_metadata <- read.csv("genome_metadata.csv")
sample_table_data <- read.csv("sample_info_table.csv") %>% 
  left_join(genome_metadata %>% select(abbrev, gen_time, length_m, lifespan) %>% group_by(abbrev)
            %>% slice_head() %>% ungroup,
            by = c("abbreviation" = "abbrev"))
reference_table_data <- read.csv("reference_info_table.csv")
cet_het <- read.csv("cetacean_genome_heterozygosity.csv") %>% 
  dplyr::select(Species, Observed_pi)

#### Load whole genome heterozygosity, homozygosity, and alternate alleles -----
data_files <- read.delim("cetacean_snps_au23.txt", sep = " ", header = FALSE) %>% 
  dplyr::rename("abbrev" = 1, "reference" = 2, "total" = 3, "hets" = 4, "hom" = 5, "aallele" = 6) %>%
  left_join(genome_metadata, by = c("abbrev", "reference")) %>% 
  left_join(cet_het, by = c("species_latin" = "Species"))

#### Summarize heterozygosity across cetaceans ---------------------------------

het_summary <- cet_het %>% 
  filter(!is.na(Observed_pi)) %>% 
  #filter(Species != "Kogia_breviceps") %>% 
  summarise(mean_pi = mean(Observed_pi), max_pi = max(Observed_pi), min_pi = min(Observed_pi)) 
  
  
#### Calculate Mutation Rate----------------------------------------------------

# data frame
whole_mrate_data <- data.frame()

for (i in 1:nrow(data_files)) {
  AALLELE=data_files$aallele[i]
  TOTAL=data_files$total[i] - AALLELE
  HET=data_files$hets[i]
  HOMALT=data_files$hom[i]
  HOMREF=TOTAL-(HOMALT+HET)
  
  # Heterozygosity
  PI=HET/TOTAL
  
  # Divergence (number of alt alleles in hets and hom. alt divided by the total number of 
  # alleles)
  DIV=(HET+(2*HOMALT))/(2*TOTAL)
  
  #Generations to coalescence
  T_TT=data_files$div_timetree[i]/data_files$gen_time[i]
  T_McG=data_files$div_mcgowen2020[i]/data_files$gen_time[i]
  T_LL=data_files$div_lloyd2021[i]/data_files$gen_time[i]
  
  #noPi mutation rate per generation
  noPI_TTgen <- DIV/(2*T_TT)
  noPI_McGgen <- DIV/(2*T_McG)
  noPI_LLgen <- DIV/(2*T_LL)
  
  #per-species Pi mutation rate per generation
  PI_TTgen <- (DIV-PI)/(2*T_TT)
  PI_McGgen <- (DIV-PI)/(2*T_McG)
  PI_LLgen <- (DIV-PI)/(2*T_LL)
  
  #noPi mutation rate per year
  noPI_TTyear <- DIV/(2*data_files$div_timetree[i])
  noPI_McGyear <- DIV/(2*data_files$div_mcgowen2020[i])
  noPI_LLyear <- DIV/(2*data_files$div_lloyd2021[i])
  
  #per-species Pi mutation rate per year
  PI_TTyear <- (DIV-PI)/(2*data_files$div_timetree[i])
  PI_McGyear <- (DIV-PI)/(2*data_files$div_mcgowen2020[i])
  PI_LLyear <- (DIV-PI)/(2*data_files$div_lloyd2021[i])
  
  #PI sensitivity - mutation rate per year with min, max, mean pi using McG divergence
  PImin_McGyear <- (DIV-het_summary$min_pi)/(2*data_files$div_mcgowen2020[i])
  PImax_McGyear <- (DIV-het_summary$max_pi)/(2*data_files$div_mcgowen2020[i])
  PImean_McGyear <- (DIV-het_summary$mean_pi)/(2*data_files$div_mcgowen2020[i])
  
  #PI sensitivity - mutation rate per generation with min, max, mean pi using McG divergence
  PImin_McGgen <- (DIV-het_summary$min_pi)/(2*T_McG)
  PImax_McGgen <- (DIV-het_summary$max_pi)/(2*T_McG)
  PImean_McGgen <- (DIV-het_summary$mean_pi)/(2*T_McG)
  
  temp_wga_data <- data.frame(species = data_files$species[i],
                              abbrev = data_files$abbrev[i],
                              Pimethod = c(rep("noPI",6),rep("PI",6), 
                                           rep(c("PIsummin","PIsummax", "PIsummean"),2)), 
                              time_step = c(rep("gen",3),rep("year",3),
                                            rep("gen",3),rep("year",3),
                                            rep("gen",3),rep("year",3)),
                              div_est = c(rep(c("timetree","McG","Lloyd"),4),rep("McG",6)),
                              rate = c(noPI_TTgen, noPI_McGgen, noPI_LLgen,
                                       noPI_TTyear, noPI_McGyear, noPI_LLyear,
                                       PI_TTgen, PI_McGgen, PI_LLgen,
                                       PI_TTyear, PI_McGyear, PI_LLyear,
                                       PImin_McGgen, PImax_McGgen, PImean_McGgen,
                                       PImin_McGyear, PImax_McGyear, PImean_McGyear), 
                              reference = data_files$reference[i],
                              pi = PI)
  
  whole_mrate_data <- rbind(whole_mrate_data, temp_wga_data)
}

#### Test for effect of pi -----------------------------------------------------

pi_test_data <- whole_mrate_data %>% 
  filter(div_est == "McG") %>% 
  filter(reference == "Pmac") %>% 
  filter(time_step == "year")

pi_normality <- shapiro.test(pi_test_data$rate)
pi_normality

pi_kw <- kruskal.test(rate ~ Pimethod, data = pi_test_data)
pi_kw

pi_test_nomax <- pi_test_data %>% 
  filter(Pimethod != "PIsummax")

# post doc dunn test
pi_dunn <- data.frame(dunn.test(pi_test_data$rate, pi_test_data$Pimethod,
                                   method = "bonferroni")) %>% 
  filter(P.adjusted < 0.05) %>% 
  separate(comparisons, into = c("pi1", "pi2"), sep = " - ")

pi_dunn_diff1 <- pi_dunn %>% 
  group_by(pi1) %>% 
  summarise(n_sig1 = n())

pi_dunn_diff <- pi_dunn %>% 
  group_by(pi2) %>% 
  summarise(n_sig2 = n()) %>% 
  full_join(pi_dunn_diff1, by = c("pi2" = "pi1")) %>% 
  mutate(across(c(n_sig1,n_sig2), 
                ~ case_when(is.na(.) ~ 0, TRUE ~ .))) %>% 
  mutate(tot_sig = n_sig1 + n_sig2) %>% 
  filter(tot_sig > 2) %>% 
  mutate(y = 2.1e-10)

ggplot(data = pi_test_data) +
  geom_violin(aes(x=Pimethod, y=rate, fill = Pimethod, color = Pimethod), alpha = 0.6) +
  geom_point(pi_dunn_diff, mapping = aes(x=pi2, y = y), shape = 8, size = 1.5, stroke = 1.1, color = "indianred3") +
  theme_light() +
  labs(x= expression(paste(pi[anc], " correction method")), y="Mutations/site/year") +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=pnw_palette(n=5, name="Shuksan")) +
  scale_color_manual(values=pnw_palette(n=5, name="Shuksan")) +
  scale_alpha(guide="none") +
  geom_dotplot(aes(x=Pimethod, y=rate), binaxis = 'y', fill = "black", color ="black", alpha = 0.5, stackdir = "center") +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0)), legend.position = "none") +
  scale_x_discrete(labels = c("noPI" = expression(paste("no ",pi)), 
                              "PI" = expression(paste("species ",pi)), 
                              "PIsummax" = expression(paste("max ",pi)), 
                              "PIsummean" = expression(paste("mean ", pi)), 
                              "PIsummin" = expression(paste("min ", pi))))
  
#### Test for effect of divergence date ----------------------------------------

div_test_data <- whole_mrate_data %>% 
  filter(reference == "Pmac") %>% 
  filter(time_step == "year") %>% 
  filter(Pimethod == "noPI")

div_normality <- shapiro.test(div_test_data$rate)
div_normality

div_anova <- aov(rate ~ div_est, data = div_test_data)
summary(div_anova)

ggplot(data = div_test_data, aes(x=div_est, y=rate, fill = div_est, color = div_est, alpha = 0.6)) +
  geom_hline(aes(yintercept=2.2e-10), linetype= "dotdash") +
  geom_violin() +
  theme_light() +
  xlab("Divergence date estimate") +
  ylab("Mutations/site/year") +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=pnw_palette(n=5, name="Winter")) +
  scale_color_manual(values=pnw_palette(n=5, name="Winter")) +
  scale_alpha(guide="none") +
  geom_dotplot(binaxis = 'y', fill = "black", color ="black", alpha = 0.5, stackdir = "center") +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0)), legend.position = "none")


###DATA FRAMES### --------------------------------------------------------------
#mutation rate by species per generation 
species_mrate_gen <- whole_mrate_data %>% 
  left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  filter(reference == "Pmac") %>% 
  filter(div_est == "McG") %>% 
  filter(time_step == "gen") %>% 
  filter(Pimethod == "noPI")

#mutation rate by species per year 
species_mrate_year <- whole_mrate_data %>% 
  left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  filter(reference == "Pmac") %>% 
  filter(div_est == "McG") %>% 
  filter(time_step == "year") %>% 
  filter(Pimethod == "noPI")

#mutation rate by species per generation w multiple references 
mult_refs_gen <- whole_mrate_data %>% 
  left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  group_by(species) %>% 
  filter(n()>18) %>% ungroup() %>%
  filter(div_est == "McG") %>%
  filter(time_step == "gen") %>% 
  filter(Pimethod == "noPI") %>% 
  mutate(ref_dist = case_when(reference == "Pmac" ~ "Intermediate",
                              reference %in% c("Egla","Oorc") ~ "Recent",
                              reference == "Hamb" ~ "Deep",
                              TRUE ~ NA))

#mutation rate by species per year w multiple references 
mult_refs_year <- whole_mrate_data %>% 
  left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  group_by(species) %>% 
  filter(n()>18) %>% ungroup() %>%
  filter(div_est == "McG") %>%
  filter(time_step == "year") %>% 
  filter(Pimethod == "noPI") %>% 
  mutate(ref_dist = case_when(reference == "Pmac" ~ "Intermediate",
                              reference %in% c("Egla","Oorc") ~ "Recent",
                              reference == "Hamb" ~ "Deep",
                              TRUE ~ NA))

#mutation rate by species w multiple references 
mult_refs <- whole_mrate_data %>% 
  left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  group_by(species) %>% 
  filter(n()>18) %>% ungroup() %>%
  filter(div_est == "McG") %>%
  filter(Pimethod == "noPI") %>% 
  mutate(ref_dist = case_when(reference == "Pmac" ~ "Intermediate",
                              reference %in% c("Egla","Oorc") ~ "Recent",
                              reference == "Hamb" ~ "Deep",
                              TRUE ~ NA))

#test whether mutation rate is normally distributed
# all species, one reference, mutation/year
hist(species_mrate_year$rate)
ggqqplot(species_mrate_year$rate)

species_rate_normality <- shapiro.test(species_mrate_year$rate)
species_rate_normality

# multi-ref species, mutation/gen
hist(mult_refs_gen$rate)
ggqqplot(mult_refs_gen$rate)

ref_normality <- shapiro.test(mult_refs_gen$rate)
ref_normality

#is mutation rate significantly different among infraorder or family?
infraorder_anova_gen <- aov(rate ~ infraorder, data = species_mrate_gen)
summary(infraorder_anova_gen)

infraorder_anova_gen_sub <- aov(rate ~ infraorder, data = species_mrate_gen %>% 
                                   filter(abbrev != "Igeo"))
summary(infraorder_anova_gen_sub)

infraorder_anova_year <- aov(rate ~ infraorder, data = species_mrate_year)
summary(infraorder_anova_year)

infraorder_anova_year_sub <- aov(rate ~ infraorder, data = species_mrate_year %>% 
                               filter(abbrev != "Igeo"))
summary(infraorder_anova_year_sub)

family_data <- species_mrate_year %>% 
  group_by(family) %>% 
  filter(n() > 1)

family_anova <- aov(rate ~ family, data = family_data)
summary(family_anova)

family_data_gen <- species_mrate_gen %>% 
  group_by(family) %>% 
  filter(n() > 1)

family_anova_gen <- aov(rate ~ family, data = family_data_gen)
summary(family_anova_gen)

#what is the range in mutation rate variability among species?

species_rate_range <- species_mrate_year %>% 
  summarise(ratediff = max(rate) - min(rate))

#is mutation rate significantly different among references?
reference_anova <- aov(rate ~ ref_dist, data = mult_refs_year)
summary(reference_anova)

reference_anova_sp <- aov(rate ~ reference, data = mult_refs_year)
summary(reference_anova_sp)

reference_anova_near <- aov(rate ~ reference, data = mult_refs_year %>% filter(ref_dist=="Recent"))
summary(reference_anova_near)

reference_anova_far <- aov(rate ~ infraorder, data = mult_refs_year %>% filter(ref_dist=="Deep"))
summary(reference_anova_far)

###PLOTS####--------------------------------------------------------------------
#plot for whole mutation rate for only sperm whale refs-------------------------
ggplot(data=species_mrate_gen, aes(x=abbrev,y=rate,shape=div_est,color=species)) +
  geom_point(size=3.7) +
  labs(x="Species", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate", 
       subtitle = "Aligned to a sperm whale reference genome with the McGowan et al (2020)\ndivergence estimates") +
  theme_light() +
  theme(text = element_text(size = 20), plot.subtitle = element_text(size = 14)) +
  guides(color="none") +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_shape_discrete(guide="none") +
  scale_color_manual(values=pnw_palette(n=14,name="Sunset2")) +
  geom_hline(aes(yintercept=1.08e-8, linetype="B")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

#plot for mutation rates of multiple references---------------------------------
ggplot(data=mult_refs_gen, aes(x=species,y=rate,shape=div_est,color=ref_dist)) +
  geom_point(size=3.7) +
  labs(x="Species", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate",
       subtitle = "Species aligned to alternate reference genomes based on phylogenetic distance\nand alternate divergence estimates") +
  theme_light() +
  theme(text = element_text(size = 20), plot.subtitle = element_text(size = 14)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_shape_discrete(name="Divergence Date Estimate") +
  scale_color_manual(values=pnw_palette(n=3,name="Sunset2"), 
                     name="Reference Distance") +
  geom_hline(aes(yintercept=1.08e-8, linetype="B")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

#plot for rate by lifespan (only Pmac refs)-------------------------------------
#linear plot
ggplot(data=species_mrate_year, aes(x=lifespan,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Lifespan", y="Mutations/site/year") +
  theme_light() +
  theme(text = element_text(size = 16)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

#plot for rate by generation time (only Pmac refs)-------------------------------------
#linear plot
ggplot(data=species_mrate_year, aes(x=gen_time,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Generation time", y="Mutations/site/year") +
  theme_light() +
  theme(text = element_text(size = 16)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

#violin of mysticetes v. odontocete whole mutation rates (only Pmac refs)-------

ggplot(data=species_mrate_year, aes(x=infraorder, y=rate)) +
  geom_violin() +
  labs(x="Infraorder", y="Mutations/site/year",
       title="Whole Genome Mutation Rate by Infraorder") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_fill_manual(values=pnw_palette(n=2, name="Sunset2")) +
  geom_point(size=2.7, alpha=0.75)

#violin of whole mutation rates by family (only Pmac refs)----------------------

ggplot(data=family_data, aes(x=family, y=rate)) +
  geom_violin(aes(fill=family), color = NA, linewidth=0, alpha = 0.5) +
  labs(x="Family", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Infraorder") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_fill_manual(values=pnw_palette(n=3, name="Sunset2")) +
  geom_point(size=3, alpha=0.9,aes(color=species))

#violin of different reference genome whole mutation rates---------------------

ggplot(data=mult_refs_gen, aes(x=reference, y=rate)) +
  geom_violin() +
  labs(x="Reference", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Reference") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  geom_point(size=2.7, alpha=0.75)

#save dataframes
save(species_mrate_gen, species_mrate_year, whole_mrate_data, mult_refs_gen, 
     sample_table_data, reference_table_data, species_rate_normality, ref_normality, reference_anova,
     family_anova, reference_anova_near, reference_anova_sp, reference_anova_far, 
     mult_refs, mult_refs_year, div_anova, div_test_data, pi_kw, pi_dunn_diff, pi_test_data, family_data, family_data_gen,
     infraorder_anova_year, infraorder_anova_gen, infraorder_anova_gen_sub, infraorder_anova_year_sub,
     family_anova_gen, het_summary, pi_normality, div_normality,species_rate_range, pi_test_nomax,
     file = "whole_mutation_rate_df.Rdata")

