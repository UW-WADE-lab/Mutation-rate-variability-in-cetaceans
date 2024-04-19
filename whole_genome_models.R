#whole genome mutation rate linear models testing and graphs
#Sophie Garrote 4/18/2024

--------------------------------------------------------------------------------
library(tidyverse)
library(gam)

setwd("C:/Users/Intern/Downloads")

#load data from whole_genome_rate R script
load("Mutation-rate-variability-in-cetaceans/whole_mutation_rate_df.Rdata")
genome_metadata <- read.csv("R_script_info_mrate.csv")

##CORRELATION TESTING-----------------------------------------------------------
# rate and lifespan
cor(whole_mrate_graph$rate, whole_mrate_graph$lifespan, use = "complete.obs")
## 0.6834166

# rate and body mass
cor(whole_mrate_graph$rate, whole_mrate_graph$body_mass, use = "complete.obs")
## 0.4883015

# rate and gen time
cor(whole_mrate_graph$rate, whole_mrate_graph$gen_time, use = "complete.obs")
## 0.9839686

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

