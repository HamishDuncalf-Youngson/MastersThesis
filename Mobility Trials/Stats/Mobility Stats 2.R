rm(list=ls()) 


library(ggplot2)
library(tidyverse)
library(emmeans)
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


##Set explanatory variables as foactors 
M$Temperature <- factor(M$Temperature)
M$Chemical <- factor(M$Chemical)
M$Chem.conc <- factor(M$Chem.conc)
M$Concentration <- factor(M$Concentration)





###############################Three Way Anova ############

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


pwcchem <- M %>%
  group_by (Temperature) %>%
  emmeans_test(Total.Distance ~ Chemical, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) 

pwcchem


write_csv(pwcchem, "/Users/hamishyoungson/Documents/Imperial/Thesis/Thesis files/Data/Mobility Trials/Stats/pwcchem.csv")



####


