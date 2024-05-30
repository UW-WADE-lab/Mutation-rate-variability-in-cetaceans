#whole genome mutation rate linear models testing and graphs
#Sophie Garrote 4/18/2024

--------------------------------------------------------------------------------
install.packages("lme4", type = "source")
install.packages("lmerTest")
install.packages("modelsummary")
install.packages("corrplot")
install.packages("ltm")
library(ltm)
library(modelsummary)
library(lme4)
library(lmerTest)
library(tidyverse)
library(gam)
library(ggplot2)
library(PNWColors)
library(corrplot)

setwd("C:/Users/Intern/Downloads")

#load data from whole_genome_rate R script
load("Mutation-rate-variability-in-cetaceans/whole_mutation_rate_df.Rdata")
genome_metadata <- read.csv("R_script_info_mrate.csv")

##CORRELATION TESTING-----------------------------------------------------------
# correlation matrix
mrate_matrix_variables <- whole_mrate_graph %>% select(c("rate","gen_time","body_mass","lifespan"))
mrate_var <- cor(mrate_matrix_variables, use = "complete.obs")

##         rate      gen_time  body_mass lifespan
#rate      1.0000000 0.9804021 0.4883015 0.7022131
#gen_time  0.9804021 1.0000000 0.5320043 0.6950965
#body_mass 0.4883015 0.5320043 1.0000000 0.8390742
#lifespan  0.7022131 0.6950965 0.8390742 1.0000000

corrplot(mrate_var, type = "upper", order = "hclust", addCoef.col = "black")

# correlation between infraorder and lifespan?
biserial.cor(whole_mrate_graph$lifespan, whole_mrate_graph$infraorder)
# 0.5760454

# Is mrate NORMALLY DISTRIBUTED?------------------------------------------------

# mrate histogram (only Pmac + McGowan)
hist(whole_mrate_graph$rate)
shapiro.test(whole_mrate_graph$rate)
# p-value = 0.63

# mrate histogram (all references + estimates)
hist(whole_mrate_data$rate)
shapiro.test(whole_mrate_data$rate)
# p-value = 1.789e-05

## Yes, mrate is normally distributed for the set of data being used in 
## modelling (only Pmac reference + McGowan estimate/whole_mrate_graph data frame)

##LIFESPAN + RATE MODELS--------------------------------------------------------

#LINEAR model and plot
lifespan_model_lm <- lm(formula = rate ~ lifespan, data=whole_mrate_graph)
summary(lifespan_model_lm)

lm.1 <- lm(rate ~ lifespan, whole_mrate_graph)
plot(lm.1$fitted, lm.1$residuals)

ggplot(data=whole_mrate_graph, aes(x=lifespan,y=rate)) +
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

#NONLINEAR model
gam.1 <- gam(rate ~ lifespan, data=whole_mrate_graph)
gam.2 <- gam(rate ~ lo(lifespan), data=whole_mrate_graph)

summary(gam.1)
summary(gam.2)

plot(gam.1$fitted, gam.1$residuals)
plot(gam.2$fitted, gam.1$residuals)

#LOESS/CURVED model and plot
lifespan_model_loess <- loess(formula = rate ~ lifespan, whole_mrate_graph, span=0.8)
summary(lifespan_model_loess)

ggplot(data=whole_mrate_graph, aes(x=lifespan,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Lifespan (years)", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Lifespan") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method=loess, color = "black", size = 0.75, alpha = 0.20)

#EXPONENTIAL model
lifespan_model_exponential <- lm(formula = log(rate) ~ lifespan, data = whole_mrate_graph)
summary(lifespan_model_exponential)

#LINEAR model removing and predicting Mden lifespan
whole_mrate_noMden <- whole_mrate_graph %>% filter(abbrev!='Mden')

lifespan_noMden <- lm(formula = rate ~ lifespan, data=whole_mrate_noMden)
summary(lifespan_noMden)

Mden_mrate <- data.frame(rate = c(1.698803e-09))

cc <- coef(lifespan_noMden)
xnew <- (((1.698803e-09)-cc[1])/cc[2])
xnew

predict(lifespan_noMden, newdata=Mden_mrate)

## Mden lifespan = 26.64755 (when function is inverted)

##LINEAR MIXED MODEL------------------------------------------------------------
# w/ a single random effect for INFRAORDER
lme.1 <- lmer(rate ~ lifespan + (1 |infraorder), whole_mrate_graph)
summary(lme.1)

lme.2 <- lm(rate ~ lifespan, whole_mrate_graph)
summary(lme.2)

AIC_lmm <- AIC(lme.1,lme.2)

modelplot(lme.1)

plot(coef(lme.1))

##BODY MASS + RATE MODELS-------------------------------------------------------
#LINEAR model
bodymass_model_lm <- lm(formula = rate ~ body_mass, data=whole_mrate_graph)
summary(bodymass_model_lm)

lm.2 <- lm(rate ~ body_mass, whole_mrate_graph)
plot(lm.2$fitted, lm.2$residuals)

ggplot(data=whole_mrate_graph, aes(x=body_mass,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Body Mass (kg)", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Body Mass") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

#LINEAR model w/ Bmus outlier removed
whole_mrate_noBmus <- whole_mrate_graph %>% filter(body_mass < 135000.000)

bodymass_noBmus <- lm(formula = rate ~ body_mass, data=whole_mrate_noBmus)
summary(bodymass_noBmus)

ggplot(data=whole_mrate_noBmus, aes(x=body_mass,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Body Mass (kg)", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Body Mass") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)
 
##GEN TIME + RATE MODELS--------------------------------------------------------
#LINEAR model
gentime_model_lm <- lm(formula = rate ~ gen_time, data=whole_mrate_graph)
summary(gentime_model_lm)

lm.3 <- lm(rate ~ gen_time, whole_mrate_graph)
plot(lm.3$fitted, lm.3$residuals)

ggplot(data=whole_mrate_graph, aes(x=gen_time,y=rate)) +
  geom_point(aes(color=factor(infraorder)), size=3.7) +
  labs(x="Generation Time (year)", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate by Generation Time") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  scale_color_manual(values=pnw_palette(n=2,name="Sunset2"),
                     name = "Infraorder") +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

##ODONTOCETE LINEAR LIFESPAN----------------------------------------------------
# dataframe for odontocetes
odontocete_mrate <- whole_mrate_graph %>%  
  filter(infraorder == "odontocete")

#linear model
odontocete_lm <- lm(formula = rate ~ lifespan, data=odontocete_mrate)
summary(odontocete_lm)
plot(odontocete_lm$fitted, odontocete_lm$residuals)

ggplot(data=odontocete_mrate, aes(x=lifespan,y=rate)) +
  geom_point( size=3.7) +
  labs(x="Lifespan", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate of Odontocetes by Lifespan") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

##MYSTICETE LINEAR LIFESPAN-----------------------------------------------------
# dataframe for mysticetes
mysticete_mrate <- whole_mrate_graph %>%  
  filter(infraorder == "mysticete")

#linear model
mysticete_lm <- lm(formula = rate ~ lifespan, data=mysticete_mrate)
summary(mysticete_lm)
plot(mysticete_lm$fitted, mysticete_lm$residuals)

ggplot(data=mysticete_mrate, aes(x=lifespan,y=rate)) +
  geom_point( size=3.7) +
  labs(x="Lifespan", y="Mutations/site/generation",
       title="Whole Genome Mutation Rate of Odontocetes by Lifespan") +
  theme_light() +
  theme(text = element_text(size = 20)) +
  theme(plot.margin = unit(c(10,30,0,0), 'pt'), axis.title.y = element_text(margin = margin(t=0,r=12,b=0,l=5)),
        axis.title.x = element_text(margin = margin(t=12,r=0,b=5,l=0))) +
  geom_smooth(method = "lm", color = "black", size = 0.75, alpha = 0.20)

##save dataframes---------------------------------------------------------------
save(lifespan_model_lm, lme.1, lme.2, AIC_lmm, file = "C:/Users/Intern/Downloads/Mutation-rate-variability-in-cetaceans/mrate_models_df.Rdata")
