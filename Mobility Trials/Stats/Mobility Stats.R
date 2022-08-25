rm(list=ls()) 

install.packages("AICcmodavg")
install.packages("stargazer")
install.packages('htmlTable')
install.packages('sjPlot')

library(tidyverse)
library(lme4)
library(stargazer)
library(htmlTable)
library(sjPlot)
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(kableExtra)
library(formattable)
library(emmeans)
library(tidyverse)
library(lme4)
library(stargazer)
library(htmlTable)
library(sjPlot)
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(kableExtra)
library(formattable)
library(rstatix)
library(report)


#read in file
M <- read.csv('stats.csv') 



#Scatter of Temperature and Chemical
par(mfrow=c(1,1))
plot(Total.Distance ~ Temperature,data=M,pch=19,cex=1.5,
     xlab=list("Temperature",cex=1.2),
     ylab = list("Total.Distance"))
abline(mobility.mod)


##Set explanatory variables as foactors 
M$Temperature <- factor(M$Temperature)
M$Chemical <- factor(M$Chemical)
M$Chem.conc <- factor(M$Chem.conc)
M$Concentration <- factor(M$Concentration)

str(M)
#Linear Model
mobility.mod <- (lm(Total.Distance ~ Temperature*Chemical*Concentration, data= M))
names(mobility.mod) 

#Residual Diagnostic Plots for Linear Model
par(mfrow=c(2,3))  #1. Linear relationship shown (homoskedastic), 2. Normally distributed (6 might be an issue?) 
##3. There is equal variance 4. No influential outliers
plot(mobility.mod)


#histogram of residuals 
hist(mobility.mod$residuals)

#Anova Table listing 
anova(mobility.mod)

#The effect of temperature interacts with chemcials (synergistic)..

#Summary table
summary(mobility.mod)



#####Two-Way Anova 
two.way <- aov(Total.Distance ~ Chemical + Temperature, data = M)
interaction <- aov(Total.Distance ~ Chemical * Temperature, data = M)


#Model Comparison

library(AICcmodavg)

model.set <- list(two.way, interaction)
model.names <- c("two.way", "interaction")

aictab(model.set, modnames = model.names)


##Interaction is model of best fit 
##Summary 
summary(interaction)


#Post-hoc test
TukeyHSD(interaction,which="Temperature")
TukeyHSD(interaction,which="TemperatureChemical")




#perform three-way ANOVA
threewayanova <- aov(Total.Distance ~ Temperature * Chemical * Concentration, data=M)

plot(threewayanova)

#view summary of three-way ANOVA
summary(threewayanova)

#Posthoc
TukeyHSD(threewayanova)


#####Linear Mixed Model 
hist(M$Total.Distance)
par(mfrow=c(1,1))

#Standardize the explanatory variables 
M$Temperature2 <- scale(M$Temperature, center = TRUE, scale = TRUE)
M$Concentration2 <- scale(M$Concentration, center = TRUE, scale = TRUE)

##Linear Model Temperature
lmtemp<- lm(Total.Distance ~ Temperature2, data = M)
summary(lmtemp)
 
#Plot
(prelim_plot <- ggplot(M, aes(x = Temperature2, y = Total.Distance)) +
    geom_point() +
    geom_smooth(method = "lm"))


##Linear Model Concentration  
lmconc <- lm(Total.Distance ~ Concentration2,data=M)
summary(lmconc)


#plot 
(prelim_plot <- ggplot(M, aes(x = Concentration2, y = Total.Distance)) +
    geom_point() +
    geom_smooth(method = "lm"))


#residuals
plot(lmtemp, which = 1) 
plot(lmconc,which=1)

plot(lmtemp, which = 2) 
plot(lmconc,which=2)

#Boxplots
boxplot(Total.Distance ~ Chemical, data = M)
boxplot(Total.Distance~ Temperature, data=M)
boxplot(Total.Distance ~ Concentration, data=M)

##Colourplot
(colour_plot <- ggplot(M, aes(x = Chemical, y = Total.Distance, colour = Temperature)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none"))


#Facet plot
(split_plot <- ggplot(aes(Chemical, Total.Distance, colour=Concentration), data = M) + 
    geom_point() + 
    facet_wrap(~ Temperature) +
    xlab("Chemical") + 
    ylab("Total Distance"))


#Modify the model 
Lm2 <- lm(Total.Distance ~ Temperature2 + Chemical, data = M)
summary(Lm2)


##Adding multiple fixed models and random effects 
mixed.lmer <- lmer(Total.Distance ~ Temperature2*Chemical*Concentration+ (1|Replicate), data = M)
summary(mixed.lmer)
anova(mixed.lmer)

## Calculating variance explained by random effects (replicate) = ~2.88% of the data explained by replicate
3.511/(3.511 + 118.220) 

#Residuals
plot(mixed.lmer)
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))

star <- stargazer(mixed.lmer, type = "html",style = "qje",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

htmlTable(star)

tab_model(mixed.lmer)



##Fipronil 
Fipronil <- subset(M, Chemical == 'Fipronil') 

##Linear Mixed Model
Fiplmer <- lmer(Total.Distance ~ Temperature*Concentration + (1|Replicate), data= Fipronil)
summary(Fiplmer)
plot(Fiplmer)

##Linear Model
Fipmod <- lm(Total.Distance ~ Temperature*Concentration, data= Fipronil)
anova(Fipmod)
plot(Fipmod)

par(mfrow=c(2,3))  
plot(FipMod)
hist(FipMod$residuals)

anova(Fipmod)

tab_model(FipMod)







#Linear Model
lm3 <- lm(Total.Distance ~ Temperature*Chem.conc, data= M)


#Residual Diagnostic Plots for Linear Model
par(mfrow=c(2,3))  #1. Linear relationship shown (homoskedastic), 2. Normally distributed (6 might be an issue?) 
##3. There is equal variance 4. No influential outliers
plot(lm3)


#histogram of residuals 
hist(lm3$residuals)

#Anova Table listing 
anova(lm3)

#The effect of temperature interacts with chemcials (synergistic)..

#Summary table
summary(lm3)


######Two-Way Anova 
two.way <- aov(Total.Distance ~ Chem.conc + Temperature, data = M)
interaction <- aov(Total.Distance ~ Chem.conc * Temperature, data = M)



#Model Comparison



model.set <- list(two.way, interaction)
model.names <- c("two.way", "interaction")

aictab(model.set, modnames = model.names)


##Interaction is model of best fit 
##Summary 

summary(interaction)


#Post-hoc test
TukeyHSD(interaction,which='Chem.conc')
TukeyHSD(interaction,which='Temperature')



###############################Two Way Anova ############

ggboxplot(M, x = "Chemical", y = "Total.Distance", 
          color = "Temperature", palette = c("red", "black","Blue", "Green"), facet.by = "Concentration")


#Diagnostics- Normality 
modelmob  <- lm(Total.Distance ~ Temperature*Chemical*Concentration, data = M)

ggqqplot(residuals(modelmob))




#Shapiro-Wilks 
shapiro_test(residuals(modelmob))

#Not significant so we can assume normalit



par(mfrow=c(1,1))
plot(lm(Total.Distance~Temperature*Chemical*Concentration, data=M), 1)

#Some points violate the homogeneity of variances, however the homogeneity has been proved in the Levenes

#Anova
res.aov <- M %>% anova_test(Total.Distance ~ Temperature*Chemical)
res.aov 

write_csv(res.aov, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/res.aov.csv")


resultsmob <- aov(Total.Distance~ Temperature*Chemical, data = M)
anova(resultsmob) %>% kbl %>% kable_material_dark()
anovamob <- resultsmob

anovamob <- anova(resultsmob)
anovamob

write_csv(anovamob, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/anova.csv")


report(resultsmob)

TukeyHSD(resultsmob, which= Chemical)

#Post Hoc testing 
posthocmob  <- lm(Total.Distance ~ Temperature*Chemical, data = M)
Chemical.effect <- M %>%
  group_by(Temperature) %>%
  anova_test(Total.Distance ~ Chemical, error = posthocmob)

Chemical.effect

write_csv(thing, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/posthoc.csv")

#######
posthocmob  <- lm(Total.Distance ~ Temperature*Chemical, data = M)
Temperature.effect <- M %>%
  group_by(Chemical) %>%
  anova_test(Total.Distance ~ Temperature, error = posthocmob)

Temperature.effect

write_csv(thing, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/posthoc.csv")


#Main Effects 
temperature.effect <- M %>%
  group_by(Chemical) %>%
  anova_test(Total.Distance ~ Temperature, error = posthocmob)

temperature.effect 

write_csv(temperature.effect, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/tempeffect.csv")

#Chemical Effect
Chemical.effect <- M %>%
  group_by(Temperature) %>%
  anova_test(Total.Distance ~ Chemical, error = posthocmob)

Chemical.effect 


#Pairwise Comparisons 
M$Chemical <- factor(M$Chemical)

pwctemp <- M %>%
  group_by( Chemical) %>%
  emmeans_test(Total.Distance ~ Temperature, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) 
pwctemp

write_csv(pwc, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/pwccem.csv")


pwcchem <- M %>%
  group_by (Temperature) %>%
  emmeans_test(Total.Distance ~ Chemical, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) 

pwcchem


write_csv(pwcchem, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Stats/pwcchem.csv")

pwcconc <- M %>%
  group_by(Chemical, Temperature) %>%
  emmeans_test(Total.Distance ~ Chemical* Temperature, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) 

pwcconc


pwcconc <- M %>%
  group_by(Temperature) %>%
  emmeans_test(Total.Distance ~ Chemical* Temperature, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) 



####


