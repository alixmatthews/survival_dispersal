#### Mite Survival, Dispersal, and Days Until Death Analyses
#### Alix Matthews
#### 4 February 2022, updated 3 April 2022
#### -- I recoded the original databases (individual bird names...) and cleaned the code so that results are reproducible and easy to follow.
#### R version 3.6.3

setwd("~/Desktop/Mites/Survival/Data")

#### DAYS UNTIL DEATH ANALYSES ####
#### load libraries, versions to side
library(lme4) # 1.1-27.1
library(AICcmodavg) # 2.3-1
library(ggplot2) # 3.3.5
library(ggthemes) # 4.2.4
library(ggpubr) # 0.4.0
library(plyr) # 1.8.6
library(dplyr) # 1.0.7
library(survminer) # 0.4.9

#### read in data and view, ensure subset is correct
data <- read.csv(file = "mite_DDT_recode.csv", sep = ",")
str(data)

data<-subset(data, Bird_ID !="")
levels(data$Bird_ID)
data$Bird_ID <- factor(data$Bird_ID) # factor it again to remove fall AMREs
levels(data$Bird_ID)

#### subset Amerodectes only
data_amero<-subset(data, Mite_sp !="Proctophyllodes quadratus")
levels(data_amero$Mite_sp) # Procs still there
data_amero$Mite_sp <- factor(data_amero$Mite_sp) # factor it again to remove Procs
levels(data_amero$Mite_sp)

data_amero<-subset(data_amero, Mite_sp !="")
levels(data_amero$Mite_sp) # blanks still there
data_amero$Mite_sp <- factor(data_amero$Mite_sp) # factor it again to remove blanks
levels(data$Bird_ID)

### Days until death versus day0 mites (by mite spp and by host spp) - basic plots
ggplot(data_amero, aes(y = days_til_death, x = day0_mites, color = Mite_sp)) +
    geom_point(size=3) +
    theme_classic()

ggplot(data_amero, aes(y = days_til_death, x = day0_mites, color = Host_sp)) +
    geom_point(size=1) +
    theme_classic()

### Spearman correlation: does the # of mites on day 0 predict the days until death?
cor.test(data_amero$day0_mites, data_amero$days_til_death, method="spearman", use = "complete.obs")

#### Not really.
#### But let's view the data.
ggscatter(data_amero, x = "day0_mites", y = "days_til_death",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "# mites day 0", ylab = "Days until death")

ggscatter(data_amero, x = "day0_mites", y = "days_til_death", color = "Mite_sp",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "# mites day 0", ylab = "Days until death")

ggscatter(data_amero, x = "day0_mites", y = "days_til_death", color = "Host_sp",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "# mites day 0", ylab = "Days until death")


### Let's check data normality and for overdispersion
#### glms
#### checking normality

# density plot of day0 mites (predictor), not bell shaped
ggdensity(data_amero$day0_mites,
          xlab = "Day 0 mites")

hist(data_amero$day0_mites)

shapiro.test(data_amero$day0_mites) # W = 0.89073, p-value = 0.01159, non-normal

# qqplots, actually not too bad
ggqqplot(data_amero$day0_mites)
qqnorm(data_amero$day0_mites)

# mean vs variance. not equal, variance is much greater than the mean, so we have overdispersion, use quasipoisson
mean(data_amero$day0_mites, na.rm=TRUE) # 88.92
var(data_amero$day0_mites, na.rm=TRUE) # 5603.493


### Let's run some GLMs now with quasipoisson error distribution
#### initial null model
mod0<-glm(days_til_death ~ 1, data = data_amero)
summary(mod0)

#### day0 mites model
mod_day0 <- glm(days_til_death ~ day0_mites, family = quasipoisson, data = data_amero)
summary(mod_day0)

#### There is no effect of day0 mites on days til death

#### Is this effect due to/related to the mite spp?
mod_day0_mite_sp_add <- glm(days_til_death ~ day0_mites+Mite_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_mite_sp_add)

#### No.

#### Ok, but is there an interaction?
mod_day0_mite_sp_mul <- glm(days_til_death ~ day0_mites*Mite_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_mite_sp_mul)

#### No.

### is there an interaction with mite spp. and host spp.?
mod_day0_mite_host_sp_mul <- glm(days_til_death ~ day0_mites*Mite_sp*Host_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_mite_host_sp_mul)

#### no

### is there an interaction with mite spp. and host spp.?
mod_day0_mite_host_sp_add <- glm(days_til_death ~ day0_mites + Mite_sp + Host_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_mite_host_sp_add)

#### no

#### Let's check by host species instead
#### Is this effect due to/related to the host spp?
mod_day0_host_sp_add <- glm(days_til_death ~ day0_mites+Host_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_host_sp_add)

#### No

#### So is there an interaction?
mod_day0_host_sp_mul <- glm(days_til_death ~ day0_mites*Host_sp, family = quasipoisson, data = data_amero)
summary(mod_day0_host_sp_mul)

#### No.


### Let's do some plotting using quasipoisson - single best fit line

### first by mite species (no effect of mite species according to glms)
# can assign this DDT_mite_spp <- then save as .tiff using ggsave
# Figure used in manuscript
DDT_mite_spp <- ggplot(data_amero, aes(y = days_til_death, x = day0_mites)) +
    geom_point(aes(col = Mite_sp, shape = Host_sp), size = 3) +
    geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE, method.args = list(family = "quasipoisson")) +
    xlab("Number of Mites on Day 0") +
    ylab("Days until Death") +
    labs(color="Mite Species") +
    labs(shape = "Host Species") +
    scale_color_manual(values=c("#0072B2", "#E69F00")) +
    theme_classic() +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16, face = "italic")) +
    #theme(legend.position = "top") +
    scale_shape_discrete(labels=c("Setophaga ruticilla","Setophaga cerulea", "Protonotaria citrea")) +
    scale_y_continuous(breaks=c(5, 7, 9, 11, 13, 15, 17)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18, face = "bold"))

DDT_mite_spp

ggsave("DDT_mite_spp.tiff", DDT_mite_spp, dpi = 600, width = 8, height = 6, units = "in")


### then by host spp. (single best fit line)
ggplot(data_amero, aes(y = days_til_death, x = day0_mites)) +
    geom_point(aes(col = Host_sp), size = 3) +
    geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE, method.args = list(family = "quasipoisson")) +
    xlab("Number of Mites on Day 0") +
    ylab("Days until Death") +
    labs(color="Host Species") +
    scale_color_manual(values=c("tomato", "#0072B2", "#E69F00")) +
    theme_classic() +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16, face = "italic")) +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18, face = "bold"))


### then by host spp. (multiple best fit lines)
ggplot(data_amero, aes(y = days_til_death, x = day0_mites, color = Host_sp)) +
    geom_point() +
    geom_smooth(method = glm, alpha = 0.3, se = TRUE, method.args = list(family = "quasipoisson")) +
    xlab("Number of Mites on Day 0") +
    ylab("Days until Death") +
    labs(color="Host Species") +
    scale_color_manual(values=c("tomato", "#0072B2", "#E69F00")) +
    theme_classic() +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16, face = "italic")) +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18, face = "bold"))


### Good to move on.

#### DATA SUMMARIES ####
### By splitting the data by host, we have fewer A. ischyros for each host spp.

#### make a loop getting the average days til death by species

#### first by host species
ddt_host <- ddply(data_amero, c("Host_sp"), summarise,
                  N    = sum(!is.na(days_til_death)),
                  mean = mean(days_til_death, na.rm=TRUE),
                  sd   = sd(days_til_death, na.rm=TRUE),
                  se   = sd / sqrt(N)
)

ddt_host$lower = ddt_host$mean + ddt_host$se
ddt_host$upper = ddt_host$mean - ddt_host$se

ddt_host

#### then by mite species
ddt_mite <- ddply(data_amero, c("Mite_sp"), summarise,
                  N    = sum(!is.na(days_til_death)),
                  mean = mean(days_til_death, na.rm=TRUE),
                  sd   = sd(days_til_death, na.rm=TRUE),
                  se   = sd / sqrt(N)
)

ddt_mite$lower = ddt_mite$mean + ddt_mite$se
ddt_mite$upper = ddt_mite$mean - ddt_mite$se

ddt_mite

### estimate # of mites per host spp.
str(data_amero)
num_mites <- ddply(data_amero, c("Mite_sp"), summarise,
                   N    = sum(!is.na(day0_mites)),
                   mean = mean(day0_mites, na.rm=TRUE),
                   sd   = sd(day0_mites, na.rm=TRUE),
                   se   = sd / sqrt(N)
)

num_mites$lower = num_mites$mean + num_mites$se
num_mites$upper = num_mites$mean - num_mites$se

num_mites


#### subset Prcos only, for fun
data_proc<-subset(data, Mite_sp !="Amerodectes ischyros")
data_proc<-subset(data_proc, Mite_sp !="Amerodectes protonotaria")
data_proc<-subset(data_proc, Mite_sp !="")

levels(data_proc$Mite_sp) # Procs still there
data_proc$Mite_sp <- factor(data_proc$Mite_sp) # factor it again to remove Ameros and blanks
levels(data_proc$Mite_sp)

ddt_proc <- ddply(data_proc, c("Mite_sp"), summarise,
                  N    = sum(!is.na(days_til_death)),
                  mean = mean(days_til_death, na.rm=TRUE),
                  sd   = sd(days_til_death, na.rm=TRUE),
                  se   = sd / sqrt(N)
)

ddt_proc$lower = ddt_proc$mean + ddt_proc$se
ddt_proc$upper = ddt_proc$mean - ddt_proc$se

ddt_proc



#### SURVIVAL ANALYSES ####
#### Load libraries  # version
library(survival) # 3.2-13
library(survminer) # 0.4.9
library(plyr) # 1.8.6
library(dplyr) # 1.0.7
library(ggplot2) # 3.3.5
library(ggthemes) # 4.2.4
library(reshape2)  # 1.4.4


#### ++++ Load and adjust survival data ####
surv_data_main <- read.csv(file = "mite_survival_2019_2020_recode.csv", sep = ",")
str(surv_data_main)

# change intergers to numeric
surv_data_main$day <- as.numeric(surv_data_main$day)
surv_data_main$event <- as.numeric(surv_data_main$event)
surv_data_main$day0_mites <- as.numeric(surv_data_main$day0_mites)

surv_data<-subset(surv_data_main, bird_id !="")

levels(surv_data$bird_id)
surv_data$bird_id <- factor(surv_data$bird_id) # factor it again
levels(surv_data$bird_id)

#### ++++ Survival on entire dataset (Proc and Amero), for fun ####
surv <- Surv(time=surv_data$day, event=surv_data$event)
# surv # take a look if you want
fit_km <- survfit(surv ~ mite_species, data = surv_data)
summary(fit_km)
pairwise_survdiff(Surv(day,event) ~ mite_species, data = surv_data, p.adjust.method = "fdr", rho = 0)
surv_pvalue(fit_km)
print(fit_km, print.rmean=TRUE)
ggsurvplot(fit_km, data=surv_data, conf.int= TRUE, pval = TRUE)
ggsurvplot(fit_km, data=surv_data, conf.int= TRUE, pval = TRUE, fun = "event")
fit_cph <- coxph(surv ~ mite_species, data = surv_data)
summary(fit_cph)
ggforest(fit_cph, data = surv_data)

#### ++++ Survival of Amerodectes ####

# subset Amerodectes only
surv_data_amero<-subset(surv_data, mite_species !="Proctophyllodes quadratus")
levels(surv_data_amero$mite_species) # Procs still there
surv_data_amero$mite_species <- factor(surv_data_amero$mite_species) # factor it again to remove Procs
levels(surv_data_amero$mite_species) # Procs gone
str(surv_data_amero)

# Create survival object; event = death
surv_amero <- Surv(time=surv_data_amero$day, event=surv_data_amero$event)

# take a look if you want
# surv_amero

# Fit a Kaplan-Meier survival curve; grouping is mite species
fit_km_amero <- survfit(surv_amero ~ mite_species, data = surv_data_amero)

# Show life table; shows survival times, proportion of surviving mites at each time point, # of death events at each time point, by each group
summary(fit_km_amero)

# Print short summary of survival curve; shows # of observations, # of events, median survival with confidence limits for the median. This shows lower survival (i.e., quicker death) for A. protonotaria compared to A. ischyros
print(fit_km_amero, print.rmean=TRUE)

# Plot the survival curve as Kaplan-Meier plots. pval of a log-rank test is displayed as well! ggplot graphics. This plot has all sorts of customizations on it to show their availability.

# with Risk table
ggsurvplot(fit_km_amero, data=surv_data_amero, conf.int= TRUE, pval = TRUE, pval.method = TRUE, risk.table = TRUE, fontsize=4, risk.table.height=0.25, legend.labs=c("Amerodectes ischyros", "Amerodectes protonotaria"), legend.title = "Species", font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"), palette = c("#0072B2", "#E69F00"), xlab = "Time (days)", conf.int.alpha = 0.55)

# without Risk table
Fig_surv<-ggsurvplot(fit_km_amero, data=surv_data_amero, conf.int= TRUE, pval = FALSE, pval.method = TRUE, risk.table = FALSE, fontsize=4, legend.labs=c("Amerodectes ischyros", "Amerodectes protonotaria"), legend.title = "Species", font.x = c(18, "bold"), font.y = c(18, "bold"), font.tickslab = c(14, "plain"), palette = c("#0072B2", "#E69F00"), xlab = "Time (days)", conf.int.alpha = 0.55)

Fig_surv$plot <- Fig_surv$plot +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16, face = "italic"))

# function to use ggsave with the ggsurvplot
grid.draw.ggsurvplot <- function(x){
    survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Figure used in manuscript
# ggsave("Fig_surv_v2.tiff", Fig_surv)


# view output from log-rank test (i.e., are there differences in survival between the two groups?)
surv_pvalue(fit_km_amero) # Yes!


# calculate pairwise comparisons between group levels with corrections for multiple testing
pairwise_survdiff(Surv(day,event) ~ mite_species, data = surv_data_amero, p.adjust.method = "fdr", rho = 0)

# can use this to obtain DF and chisq value
survdiff(Surv(day,event) ~ mite_species, data = surv_data_amero, rho = 0)


# Fit a Cox PH regression model, just by mite species; this can assess both categorical and continuous variables, and can model the effect of multiple variables at once.
fit_cph_amero <- coxph(surv_amero ~ mite_species, data = surv_data_amero)

# Take a look at the output. Looking at the exp(coef) column (which is the Hazard Ratio, HR), A. protonotaria die at 1.84549x the rate per unit time as A. ischyros (A. ischyros is the baseline). HR = 1: no effect; HR>1: increase in hazard (of dying); HR<1: reduction in hazard (more protected from dying)
fit_cph_amero
summary(fit_cph_amero) # another view of the output summary

# extract the coefficients from that model
fit_cph_amero_coeffs <- coef(summary(fit_cph_amero))
fit_cph_amero_coeffs


# Illustrate a forest plot. It shows so-called hazard ratios (HR) which are derived from the model for all covariates that we included in the formula in coxph. Briefly, an HR > 1 = an increased risk of death, while HR < 1 = a decreased risk
ggforest(fit_cph_amero, data = surv_data_amero)

# This shows that A. protonotaria have a much higher risk of death than A. ischyros. A. prot HR = 1.8 (CI = 1.7-2), meaning they have an increased risk of death compared to A. ischyros, with high significance (<0.001)




#### DISPERSAL ANALYSES ####
#### Load libraries  # version
library(survival) # 3.2-13
library(survminer) # 0.4.9
library(plyr) # 1.8.6
library(dplyr) # 1.0.7
library(ggplot2) # 3.3.5
library(ggthemes) # 4.2.4
library(reshape2)  # 1.4.4

#### ++++ Load and adjust dispersal data ####
# this is the one 'with returnees'
disp_data_main <- read.csv(file = "mite_dispersal_2019_2020_recode.csv", sep = ",")
str(disp_data_main)

disp_data<-subset(disp_data_main, bird_id !="")

levels(disp_data$bird_id)
disp_data$bird_id <- factor(disp_data$bird_id) # factor it again to remove fall AMREs
levels(disp_data$bird_id)

#### ++++ Dispersal on entire dataset (Proc and Amero), for fun ####

surv <- Surv(time=disp_data$day, event=disp_data$event)
# surv, take a look if wanted
fit_km <- survfit(surv ~ mite_species, data = disp_data)
summary(fit_km)
pairwise_survdiff(Surv(day,event) ~ mite_species, data = disp_data, p.adjust.method = "fdr", rho = 0)

ggsurvplot(fit_km, data=disp_data, conf.int= TRUE, pval = TRUE)
ggsurvplot(fit_km, data=disp_data, conf.int= TRUE, pval = TRUE, fun = "event")
fit_cph <- coxph(surv ~ mite_species, data = disp_data)
summary(fit_cph)
ggforest(fit_cph, data = disp_data)



#### ++++ Dispersal of Amerodectes ####

# subset Amerodectes only
disp_data_amero<-subset(disp_data, mite_species !="Proctophyllodes quadratus")
levels(disp_data_amero$mite_species) # Procs still there
disp_data_amero$mite_species <- factor(disp_data_amero$mite_species) # factor it again to remove Procs
levels(disp_data_amero$mite_species) # Procs gone

# Create survival object (using dispersal as the event instead of death, i.e., "time-to-disperse" (or otherwise die))
surv <- Surv(time=disp_data_amero$day, event=disp_data_amero$event)

# take a look if you want
# surv

# Fit a Kaplan-Meier survival curve; grouping is mite species
fit_km <- survfit(surv ~ mite_species, data = disp_data_amero)

# Show life table; shows dispersal times, proportion of dispersing mites at each time point, # of dispersal events at each time point, by each group
summary(fit_km)

# Print short summary of dispersal curve; shows # of observations, # of events, median time-to-dispersal with confidence limits for the median. This shows lower "staying put" (i.e., quicker dispersal) for A. protonotaria compared to A. ischyros
print(fit_km)


surv_pvalue(fit_km)

pairwise_survdiff(Surv(day,event) ~ mite_species, data = disp_data_amero, p.adjust.method = "fdr", rho = 0)

# can use this to obtain DF and chisq value
survdiff(Surv(day,event) ~ mite_species, data = disp_data_amero, rho = 0)



ggsurvplot(fit_km, data=disp_data_amero, conf.int= TRUE, pval = TRUE)
ggsurvplot(fit_km, data=disp_data_amero, conf.int= TRUE, pval = TRUE, fun = "cumhaz")


Fig_disp<-ggsurvplot(fit_km, data=disp_data_amero, conf.int= TRUE, pval = FALSE, fun = "cumhaz", risk.table = FALSE, fontsize=4, legend.labs=c("Amerodectes ischyros", "Amerodectes protonotaria"), legend.title = "Species", font.x = c(18, "bold"), font.y = c(18, "bold"), font.tickslab = c(14, "plain"), palette = c("#0072B2", "#E69F00"), xlab = "Time (days)", ylab = "Cumulative dispersal hazard", conf.int.alpha = 0.55)


Fig_disp$plot <- Fig_disp$plot +
    theme(legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16, face = "italic"))

# function to use ggsave with the ggsurvplot
grid.draw.ggsurvplot <- function(x){
    survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Figure used in manuscript
# ggsave("Fig_disp2.tiff", Fig_disp)


fit_cph <- coxph(surv ~ mite_species, data = disp_data_amero)
summary(fit_cph)
ggforest(fit_cph, data = disp_data_amero)





