rm(list=ls()) 

#Dependencies
install.packages("kableExtra")
install.packages("formattable")
install.packages("rstatix")
install.packages("report")


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


#Data 
Pred <- read.csv("PredatorPrey Trials 2.csv")


#Explanatory variables as factors 
Pred$Temperature <- factor(Pred$Temperature)
Pred$Treatment <- factor(Pred$Treatment)

str(Pred)


#Two Way Anova 
#- plot 
ggboxplot(Pred, x = "Treatment", y = "Prey.eaten", color = "Temperature", palette = c("red", "black"))




#Diagnostics- Normality 
modelpred  <- lm(Prey.eaten ~ Temperature*Treatment, data = Pred)


ggqqplot(residuals(modelpred))


#Shapiro-Wilks 
shapiro_test(residuals(modelpred)) #Not significant so we can assume normality




#The group outliers shown on this graph ^

##Levenes Test
Pred %>% levene_test(Prey.eaten ~ Temperature*Treatment) 
#Levenes Not significant, so we can assume homogeneity of variances in the different groups ^


par(mfrow=c(1,1))
plot(lm(Prey.eaten~Temperature*Treatment, data=Pred), 1)
#One point violates the homogeneity of variances, however the homogeneity has been proved in the Levenes 

#Anova Test
res.aov <- Pred %>% anova_test(Prey.eaten ~ Temperature*Treatment)
res.aov

write_csv(res.aov, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Predator Prey trials/Outputs/Tables/resaov.csv")

#Anova
resultspred <- aov(Prey.eaten ~ Temperature * Treatment, data = Pred)
aresultspred<- anova(resultspred) 
anova(resultspred)

write_csv(aresultspred, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Predator Prey trials/Outputs/Tables/anova.csv")
report(resultspred)



#Post Hoc testing 
posthocpred  <- lm(Prey.eaten ~ Temperature*Treatment, data = Pred)
Pred %>%
  anova_test(Prey.eaten ~ Treatment*Temperature, error = posthocpred)


#Main Effects 
temperature.effect <- Pred %>%
  group_by(Temperature) %>%
  anova_test(Prey.eaten ~ Treatment, error = posthocpred)

temperature.effect




#Pairwise Comparisons 

pwc <- Pred %>%
  group_by(Temperature) %>%
  emmeans_test(Prey.eaten ~ Treatment, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk

pwc




write_csv(pwc, "/Users/hamishyoungson/Documents/Imperial/Thesis/Thesis files/Data/Predator Prey trials/pwc.csv")







