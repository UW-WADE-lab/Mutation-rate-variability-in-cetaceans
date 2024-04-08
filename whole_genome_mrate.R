#looping DIV and PI calculations for mutation rate calculation

#Sophie 11/21
###############################################################################

library(tidyverse)

setwd("C:/Users/Intern/Downloads")

# get genome metadata-----------------------------------------------------------

genome_metadata <- read.csv("R_script_info_mrate.csv")

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
    
    #Time to coalesence
    T_TT=data_files$div_timetree[i]/data_files$gen_time[i]
    T_McG=data_files$div_mcgowen2020[i]/data_files$gen_time[i]
    T_LL=data_files$div_lloyd2021[i]/data_files$gen_time[i]
    
    #noPi mutation rate
    noPI_TT <- DIV/(2*T_TT)
    noPI_McG <- DIV/(2*T_McG)
    noPI_LL <- DIV/(2*T_LL)
    
    #Pi mutation rate
    #PI_mutation <- (DIV-PI)/(2*T_TT)
    
    temp_wga_data <- data.frame(species = rep(data_files$species[i],3),
                                abbrev = rep(data_files$abbrev[i],3),
                                #method = c("noPI","noPI","noPI","PI"), 
                                div_est = c("timetree","McGowan","Lloyd"),
                                rate = c(noPI_TT, noPI_McG, noPI_LL), 
                                reference = rep(data_files$reference[i],3))
    
    whole_mrate_data <- rbind(whole_mrate_data, temp_wga_data)
  }
  
library(ggplot2)
library(PNWColors) #color palette library

###DATA FRAMES###---------------------------------------------------------------
# data frame for estimations of mutation rate-----------------------------------
whole_mrate_graph <- whole_mrate_data %>% left_join(data_files, by = c("abbrev", "reference")) %>% 
  filter(reference == "Pmac") %>% 
  filter(div_est == "McGowan")
est2 <- data.frame(yint=1.08e-8)

#dataframe for mutation rates of multiple references----------------------------
mult_refs_graph <- whole_mrate_data %>% group_by(species) %>% 
  filter(n()>3) %>% ungroup() %>%
  mutate(ref_dist = reference) 

mult_refs_graph$ref_dist[mult_refs_graph$ref_dist == 'Pmac'] <- 'Middle'
mult_refs_graph$ref_dist[mult_refs_graph$ref_dist == 'Egla'] <- 'Close'
mult_refs_graph$ref_dist[mult_refs_graph$ref_dist == 'Oorc'] <- 'Close'
mult_refs_graph$ref_dist[mult_refs_graph$ref_dist == 'Hamb'] <- 'Far'

est2 <- data.frame(yint=1.08e-8)

#dataframe for scatterplot with lifespan data-----------------------------------
mrate_lifespan <- whole_mrate_graph %>% 
  left_join(data_files, by = c("abbrev", "reference"))


###PLOTS####--------------------------------------------------------------------
#plot for whole mutation rate for only sperm whale refs-------------------------
ggplot(data=whole_mrate_graph, aes(x=abbrev,y=rate,shape=div_est,color=abbrev)) +
  geom_point(size=3.7) +
  labs(x="Species", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  guides(color=FALSE) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_shape_discrete(name="Divergence Date Estimate") +
  scale_color_manual(values=pnw_palette(n=12,name="Sunset2")) +
  geom_hline(data=est2, mapping=aes(yintercept=yint, linetype="B")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

#plot for mutation rates of multiple references---------------------------------
ggplot(data=mult_refs_graph, aes(x=species,y=rate,shape=div_est,color=ref_dist)) +
  geom_point(size=3.7) +
  labs(x="Species", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_shape_discrete(name="Divergence Date Estimate") +
  scale_color_manual(values=pnw_palette(n=3,name="Sunset2"), 
                     name="Reference Distance") +
  geom_hline(data=est2, mapping=aes(yintercept=yint, linetype="B")) +
  scale_linetype_discrete(labels=c('Odontocete nuclear\nmutation rate\n(Dornburg et al. 2012)'), 
                          name="Estimates")

#plot for rate by lifespan (only Pmac refs)-------------------------------------
ggplot(data=mrate_lifespan, aes(x=lifespan,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Lifespan", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Lifespan") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

#boxplot of mysticetes v. odontocete whole mutation rates (only Pmac refs)------

ggplot(data=whole_mrate_graph, aes(x=infraorder, y=rate)) +
  geom_boxplot(color="black") +
  labs(x="Infraorder", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Infraorder") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_fill_manual(values=pnw_palette(name="Sunset2"))

#boxplot of different reference genome whole mutation rates---------------------

ggplot(data=mult_refs_graph, aes(x=reference, y=rate)) +
  geom_boxplot(color="black") +
  labs(x="Reference", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Reference") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_fill_manual(values=pnw_palette(name="Sunset2"))




