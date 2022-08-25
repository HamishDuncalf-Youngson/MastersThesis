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
F <- read.csv("Fipronil concentrations.xlsx.csv")
Im <- read.csv("Imidacloprid Concentrations.csv")
P <- read.csv("P Mix All Concentrations.csv")

names(F)[2] <- "Distance"
names(F)[3] <- "Temperature"
names(F)[4] <- "Concentration"

names(Im)[2] <- "Distance"
names(Im)[3] <- "Temperature"
names(Im)[4] <- "Concentration"

names(P)[2] <- "Distance"
names(P)[3] <- "Temperature"
names(P)[4] <- "Concentration"


#Explanatory variables as factors 
F$Temperature <- factor(F$Temperature)
F$Concentration <- factor(F$Concentration)

Im$Temperature <- factor(Im$Temperature)
Im$Concentration <- factor(Im$Concentration)

P$Temperature <- factor(P$Temperature)
P$Concentration <- factor(P$Concentration)

str(F)
str(Im)
str(P)


# Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Plots


PSUM <- summarySE(P,measurevar= "Distance" , groupvars=c("Temperature", "Concentration"))

PSUM2 <- PSUM 
PSUM2$Concentration <- factor(PSUM2$Concentration)
ggplot(PSUM2, aes(x= Concentration, y=Distance, fill = factor(Temperature, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature.ºC", title = "Pesticide Mix") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se),  width=.22, position=position_dodge(.9)) + xlab("Concentration (µg/L)") +
  ylab("Distance (cm)") 





FIPSUM <- summarySE(F,measurevar= "Distance", groupvars=c("Temperature", "Concentration"))
FIPSUM

FIPSUM2 <- FIPSUM 
FIPSUM2$Concentration <- factor(FIPSUM2$Concentration)
ggplot(FIPSUM2, aes(x= Concentration, y=Distance, fill = factor(Temperature, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature.ºC", title = "Fipronil") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se),  width=.22, position=position_dodge(.9))  + xlab("Concentration (µg/L)") +
  ylab("Distance (cm)") 



ISUM <- summarySE(Im,measurevar= "Distance", groupvars=c("Temperature", "Concentration"))
ISUM

ISUM2 <- ISUM 
ISUM2$Concentration <- factor(ISUM2$Concentration)
ggplot(ISUM2, aes(x= Concentration, y=Distance, fill = factor(Temperature, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature.ºC", title = "Imidacloprid") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se),  width=.22, position=position_dodge(.9)) 


#Two Way Anova 



#Diagnostics- Normality 
modelF  <- lm(Distance ~ Temperature*Concentration, data = F)
ggqqplot(residuals(modelF))

modelIm <- lm(Distance ~ Temperature*Concentration, data = Im)
ggqqplot(residuals(modelIm))

modelP <- lm(Distance ~ Temperature*Concentration, data = P)
ggqqplot(residuals(modelP))


#Shapiro-Wilks 
shapiro_test(residuals(modelF)) #Not significant so we can assume normality
shapiro_test(residuals(modelIm)) 
shapiro_test(residuals(modelP)) 



#Anova Test
res.aovF <- F %>% anova_test(Distance ~ Temperature*Concentration)
res.aovF


write_csv(res.aovF, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/aovF.csv")

res.aovIm <- Im %>% anova_test(Distance ~ Temperature*Concentration)
res.aovIm

write_csv(res.aovIm, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/aovIm.csv")

res.aovP <- P %>% anova_test(Distance ~ Temperature*Concentration)
res.aovP

write_csv(res.aovP, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/aovP.csv")





#Main Effects 
Concentration.effectf <- F %>%
  group_by(Temperature) %>%
  anova_test(Distance ~ Concentration, error = posthocF)

Concentration.effectf

write_csv(Concentration.effectf, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/conceffectF.csv")

Concentration.effectIm <- Im %>%
  group_by(Temperature) %>%
  anova_test(Distance ~ Concentration, error = posthocF)

Concentration.effectIm

write_csv(Concentration.effectIm, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/conceffectIm.csv")


Concentration.effectP <- P %>%
  group_by(Temperature) %>%
  anova_test(Distance ~ Concentration, error = posthocF)

Concentration.effectP

write_csv(Concentration.effectP, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/conceffectP.csv")



#Pairwise Comparisons 

pwc <- P %>%
  group_by(Temperature) %>%
  emmeans_test(Distance ~ Concentration, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk

pwc

##Maybe get rid of the 'group by' above? 


pwcP <- P %>%
  group_by(Temperature) %>%
  emmeans_test(Distance ~ Concentration, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk

pwcP


write_csv(pwcP, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/pwcP.csv")



pwcF <- F %>%
  group_by(Temperature) %>%
  emmeans_test(Distance ~ Concentration, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk

pwcF


write_csv(pwcF, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/pwcF.csv")


pwcIm <- Im %>%
  group_by(Temperature) %>%
  emmeans_test(Distance ~ Concentration, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk

pwcIm


write_csv(pwcF, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Mobility Trials/Concentration Subset/pwcIm.csv")



