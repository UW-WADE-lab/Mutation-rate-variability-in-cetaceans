#whole genome mutation rate linear models testing and graphs
#Sophie Garrote 
#4/18/2024
#Modfied by AVC 7/2024

## set up environment ----------------------------------------------------------
library(ltm)
library(modelsummary)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggplot2)
library(PNWColors)
library(corrplot)

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(evomap)


#load mutation data 
load("whole_mutation_rate_df.Rdata") #from whole_genome_rate R script

#load metadata
genome_metadata <- read.csv("genome_metadata.csv")

#load cetacean phylogeny
cet_tree <- read.tree("cetacea_species.nwk")
cet_tree <- keep.tip(cet_tree, tip = genome_metadata$species_latin)

##CORRELATION TESTING-----------------------------------------------------------
# correlation matrix
demographic_variables <- species_mrate_gen %>% dplyr::select(c("gen_time","length_m","lifespan"))
demographic_correlation <- cor(demographic_variables, use = "complete.obs")
demographic_correlation

corrplot(demographic_correlation, method = "shade",type = "upper", 
         order = "hclust", addCoef.col = "black", diag = FALSE)

# correlation between infraorder and lifespan?
life_order_corr <- biserial.cor(species_mrate_gen$lifespan, species_mrate_gen$infraorder)
# 0.5760454

# Is mrate NORMALLY DISTRIBUTED?------------------------------------------------

# mrate histogram (only Pmac + McGowen)
hist(species_mrate_gen$rate)
shapiro.test(species_mrate_gen$rate)
# p-value = 0.63

# mrate histogram (all references + estimates)
hist(whole_mrate_data$rate)
shapiro.test(whole_mrate_data$rate)
# p-value = 2.2e-16

## Yes, mrate is normally distributed for the set of data being used in 
## modeling (only Pmac reference + McGowen estimate/whole_mrate_graph data frame)

##LIFESPAN + RATE MODELS--------------------------------------------------------

#Phylogenetic Generalized Least Squares (PGLS)

#mutations/year by lifespan
lifespan_pgls_year <- gls(rate ~ lifespan, 
                          correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                     data=species_mrate_year, method = "ML")

summary(lifespan_pgls_year)

pgls_year_models <- data.frame(model = c("lifespan", "length (m)", "generation time",
                                         "all parameters"),
                               AIC = c(AIC(gls(rate ~ lifespan, 
                                               correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                                               data=species_mrate_year, method = "ML")),
                                       AIC(gls(rate ~ length_m, 
                                               correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                                               data=species_mrate_year, method = "ML")),
                                       AIC(gls(rate ~ gen_time, 
                                               correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                                               data=species_mrate_year, method = "ML")),
                                       AIC(gls(rate ~ lifespan + length_m + gen_time, 
                                               correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                                               data=species_mrate_year, method = "ML"))))

#mutations/generation by lifespan
lifespan_pgls_gen <- gls(rate ~ lifespan, correlation = corBrownian(phy = cet_tree, form = ~species_latin),
                     data=species_mrate_gen, method = "REML")

summary(lifespan_pgls_gen)
AIC(lifespan_pgls_gen)

lifespan_tTable <- summary(lifespan_pgls_gen)$tTable
pGLS_ci<-gls.ci(species_mrate_gen$rate,species_mrate_gen$lifespan,vcv(cet_tree))

par(mfrow=c(1,1))
plot(rate~lifespan, data = species_mrate_gen)
abline(lifespan_pgls_gen)
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Lower2.5,lty=2)
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Upper2.5,lty=2)

ci <- data.frame(x = pGLS_ci[["CI.plot"]][["X"]], 
                 ymin = pGLS_ci[["CI.plot"]][["Lower2.5"]],
                 ymax = pGLS_ci[["CI.plot"]][["Upper2.5"]])

ggplot(data = species_mrate_gen) +
  geom_point(aes(x = lifespan, y = rate, color = infraorder)) +
  theme_light() +
  geom_abline(slope = coef(lifespan_pgls_gen)[2],
              intercept = coef(lifespan_pgls_gen)[1]) +
  geom_ribbon(data = ci, aes(x= x, ymin = ymin, 
                                  ymax = ymax), alpha = 0.4, fill = "grey50")+
  xlim(20,150)+
  ylim(0,1e-8) +
  labs(x="Lifespan", y="Mutations/site/generation") +
  theme_light() +
  theme(text = element_text(size = 16)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=6,name="Cascades")[c(3,6)],
                     name = "Infraorder")

## LINEAR model and plot -------------------------------------------------------
#lifespan
lifespan_lm <- lm(formula = rate ~ lifespan, data=species_mrate_year)
summary(lifespan_lm)

plot(lifespan_lm$fitted, lifespan_lm$residuals)

#generation time
generation_lm <- lm(formula = rate ~ gen_time, data=species_mrate_year)
summary(generation_lm)

plot(generation_lm$fitted, generation_lm$residuals)

#lifespan + generation_time
lifegen_lm <- lm(formula = rate ~ lifespan + gen_time, data = species_mrate_year)
summary(lifegen_lm)

##save dataframes---------------------------------------------------------------
save(life_order_corr, lifespan_pgls_gen, ci, lifespan_tTable, file = "mrate_models_df.Rdata")
