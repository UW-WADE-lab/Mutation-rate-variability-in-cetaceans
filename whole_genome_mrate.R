#looping DIV and PI calculations for mutation rate calculation

#Sophie 11/21
###############################################################################

library(tidyverse)
library(ggpubr)
library(PNWColors)

# get genome metadata and table data--------------------------------------------

genome_metadata <- read.csv("genome_metadata.csv")
sample_table_data <- read.csv("sample_info_table.csv")
reference_table_data <- read.csv("reference_info_table.csv")

# access files with whole genome numbers and altering column titles-------------
data_files <- read.delim("cetacean_snps_au23.txt", sep = " ", header = FALSE) %>% 
  dplyr::rename("abbrev" = 1, "reference" = 2, "total" = 3, "hets" = 4, "hom" = 5, "aallele" = 6) %>%
  left_join(genome_metadata, by = c("abbrev", "reference"))


# data frame
whole_mrate_data <- data.frame()

# LOOP FOR CALCULATING MUTATION RATE--------------------------------------------
# loop that iterates through mutation rate calculation w/ info from each row
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
  
  #Generations to coalesence
  T_TT=data_files$div_timetree[i]/data_files$gen_time[i]
  T_McG=data_files$div_mcgowen2020[i]/data_files$gen_time[i]
  T_LL=data_files$div_lloyd2021[i]/data_files$gen_time[i]
  
  #noPi mutation rate per generation
  noPI_TTgen <- DIV/(2*T_TT)
  noPI_McGgen <- DIV/(2*T_McG)
  noPI_LLgen <- DIV/(2*T_LL)
  
  #Pi mutation rate per generation
  PI_TTgen <- (DIV-PI)/(2*T_TT)
  PI_McGgen <- (DIV-PI)/(2*T_McG)
  PI_LLgen <- (DIV-PI)/(2*T_LL)
  
  #noPi mutation rate per year
  noPI_TTyear <- DIV/(2*data_files$div_timetree[i])
  noPI_McGyear <- DIV/(2*data_files$div_mcgowen2020[i])
  noPI_LLyear <- DIV/(2*data_files$div_lloyd2021[i])
  
  #Pi mutation rate per year
  PI_TTyear <- (DIV-PI)/(2*data_files$div_timetree[i])
  PI_McGyear <- (DIV-PI)/(2*data_files$div_mcgowen2020[i])
  PI_LLyear <- (DIV-PI)/(2*data_files$div_lloyd2021[i])
  
  temp_wga_data <- data.frame(species = rep(data_files$species[i],3),
                              abbrev = rep(data_files$abbrev[i],3),
                              Pimethod = c(rep("noPI",6),rep("PI",6)), 
                              time_step = c(rep("gen",3),rep("year",3),
                                            rep("gen",3),rep("year",3)),
                              div_est = rep(c("timetree","McGowan","Lloyd"),4),
                              rate = c(noPI_TTgen, noPI_McGgen, noPI_LLgen,
                                       noPI_TTyear, noPI_McGyear, noPI_LLyear,
                                       PI_TTgen, PI_McGgen, PI_LLgen,
                                       PI_TTyear, PI_McGyear, PI_LLyear), 
                              reference = rep(data_files$reference[i],3))
  
  whole_mrate_data <- rbind(whole_mrate_data, temp_wga_data)
}

###DATA FRAMES###---------------------------------------------------------------
#mutation rate by species per generation----------------------------------------
species_mrate_gen <- whole_mrate_data %>% left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  filter(reference == "Pmac") %>% 
  filter(div_est == "McGowan") %>% 
  filter(time_step == "gen") %>% 
  filter(Pimethod == "noPI")

#mutation rate by species per year----------------------------------------
species_mrate_year <- whole_mrate_data %>% left_join(data_files, by = c("abbrev", "reference", "species")) %>% 
  filter(reference == "Pmac") %>% 
  filter(div_est == "McGowan") %>% 
  filter(time_step == "year") %>% 
  filter(Pimethod == "noPI")

#mutation rate by species per generation w multiple references------------------
mult_refs_gen <- whole_mrate_data %>% group_by(species) %>% 
  filter(n()>12) %>% ungroup() %>%
  filter(div_est == "McGowan") %>%
  filter(time_step == "gen") %>% 
  filter(Pimethod == "noPI") %>% 
  mutate(ref_dist = case_when(reference == "Pmac" ~ "Middle",
                              reference %in% c("Egla","Oorc") ~ "Close",
                              reference == "Hamb" ~ "Far",
                              TRUE ~ NA))

#test whether mutation rate is normally distributed
# all species, one reference, mutation/gen
hist(species_mrate_gen$rate)
ggqqplot(species_mrate_gen$rate)

rate_normality <- shapiro.test(species_mrate_gen$rate)
rate_normality

# multi-ref species, mutation/gen
hist(mult_refs_gen$rate)
ggqqplot(mult_refs_gen$rate)

ref_normality <- shapiro.test(mult_refs_gen$rate)
ref_normality

#is mutation rate significantly different among infraorder or family?
infraorder_anova <- aov(rate ~ infraorder, data = species_mrate_gen)
summary(infraorder_anova)

family_data <- species_mrate_gen %>% 
  group_by(family) %>% 
  filter(n() > 1)

family_anova <- aov(rate ~ family, data = family_data)
summary(family_anova)

two_family_data <- family_data %>% 
  filter(family != "Delphinidae")

two_family_anova <- aov(rate ~ family, data = two_family_data)
summary(two_family_anova)

#is mutation rate significantly different among references?
reference_anova <- aov(rate ~ ref_dist, data = mult_refs_gen)
summary(reference_anova)

reference_anova_sp <- aov(rate ~ reference, data = mult_refs_gen)
summary(reference_anova_sp)

reference_anova_near <- aov(rate ~ reference, data = mult_refs_gen %>% filter(ref_dist=="Close"))
summary(reference_anova_near)

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
  labs(x="Lifespan", y="Mutations/site/year",
       title="Whole Genome Mutation Rate by Lifespan") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

#plot for rate by generation time (only Pmac refs)-------------------------------------
#linear plot
ggplot(data=species_mrate_year, aes(x=gen_time,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Lifespan", y="Mutations/site/year",
       title="Whole Genome Mutation Rate by Generation Time") +
  theme_light() +
  theme(text = element_text(size = 20)) +
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
     sample_table_data, reference_table_data, rate_normality, ref_normality, reference_anova,
     infraorder_anova, family_anova, reference_anova_near, reference_anova_sp,
     file = "whole_mutation_rate_df.Rdata")

