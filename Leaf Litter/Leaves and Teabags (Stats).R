rm(list=ls()) 

#Dependencies

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
Leaf <- read.csv("All Leaf-litter.csv")
Tea <- read.csv("All Teabags.csv")
Mort <- read.csv("Leaf litter mortality .csv")


#Calculate k values
Leaf$Leafk <- (log(Leaf$Weight/Leaf$Initial.Weight)/21)*-1
Tea$Teak <- (log(Tea$Weight..g./Tea$Initial.weight..g.)/21)*-1


#Explanatory variables as factors 
Leaf$Temperature <- factor(Leaf$Temperature)
Leaf$Treatment <- factor(Leaf$Treatment)
Leaf$Invertebrates <- factor(Leaf$Invertebrates)

Tea$Temperature <- factor(Tea$Temperature)
Tea$Treatment <- factor(Tea$Treatment)
Tea$Invertebrates <- factor(Tea$Invertebrates)

#plots 
par(mfrow=c(2,3))
plot(Leafk ~ Temperature,data=Leaf,pch=19,cex=1.5,
     xlab=list("Temperature",cex=1.2),
     ylab = list("Decay (k)"))


boxplot(Leafk ~ Treatment,data=Leaf,pch=19,cex=1.5,
     xlab=list("Chemical",cex=1.2),
     ylab = list("Decay (k)"))

boxplot(Leafk ~ Invertebrates,data=Leaf,pch=19,cex=1.5,
        xlab=list("Invertebrates",cex=1.2),
        ylab = list("Decay (k)"))


par(mfrow=c(2,3))
plot(Teak ~ Temperature,data=Tea,pch=19,cex=1.5,
     xlab=list("Temperature",cex=1.2),
     ylab = list("Decay (k)"))


boxplot(Teak ~ Treatment,data=Tea,pch=19,cex=1.5,
        xlab=list("Chemical",cex=1.2),
        ylab = list("Decay (k)"))

boxplot(Teak ~ Invertebrates,data=Tea,pch=19,cex=1.5,
        xlab=list("Invertebrates",cex=1.2),
        ylab = list("Decay (k)"))


##Linear Models

leaflm <- (lm(Leafk ~ Temperature*Treatment*Invertebrates, data= Leaf))
leaflm2 <- (lm(Leafk ~ Temperature*Treatment+Invertebrates, data= Leaf))
leaflm3 <- (lm(Leafk ~ Temperature+Treatment+Invertebrates, data= Leaf))

tealm <- (lm(Teak ~ Temperature*Treatment*Invertebrates, data= Tea))
tealm2 <- (lm(Teak ~ Temperature*Treatment+Invertebrates, data= Tea))
tealm3 <- (lm(Teak ~ Temperature+Treatment+Invertebrates, data= Tea))

#Plot Diagnostics
par(mfrow=c(2,3))
#Leaves
plot(leaflm)
hist(leaflm$residuals)

plot(leaflm2)
hist(leaflm2$residuals)

plot(leaflm3)
hist(leaflm3$residuals)

#Tea
plot(tealm)
hist(tealm$residuals)

plot(tealm2)
hist(tealm2$residuals)

plot(tealm3)
hist(tealm3$residuals)

#Model Selection 
model.set <- list(leaflm, leaflm2, leaflm3)
model.names <- c("leaflm", "leaflm2", "leaflm3")
aictab(model.set, modnames = model.names)

#leaflm model of best fit ^ 

model.set <- list(tealm, tealm2, tealm3)
model.names <- c("teaflm", "teaflm2", "teaflm3")
aictab(model.set, modnames = model.names)
#tealm model of best fit 

#Anova table 
anova(leaflm)
anova(tealm)

summary(leaflm)
summary(tealm)


#Three Way Anova - plot 
ggboxplot(Leaf, x = "Treatment", y = "Leafk", 
          color = "Temperature", palette = c("red", "black"), facet.by = "Invertebrates")

ggboxplot(Tea, x = "Treatment", y = "Teak", 
          color = "Temperature", palette = c("red", "black"), facet.by = "Invertebrates")

####Turn into barcharts? ^^^

library(kableExtra)
library(formattable)



#Diagnostics- Normality 
modelleaf  <- lm(Leafk ~ Temperature*Treatment*Invertebrates, data = Leaf)
modeltea <- lm(Teak ~ Temperature*Treatment*Invertebrates,data=Tea)

ggqqplot(residuals(modelleaf))
ggqqplot(residuals(modeltea))

#Shapiro-Wilks 
shapiro_test(residuals(modelleaf))
shapiro_test(residuals(modeltea))

#Not significant so we can assume normality

##Check normality assumption by groups.
#Leaves
sw <- Leaf %>%
  group_by(Temperature, Treatment, Invertebrates) %>%
  shapiro_test(Leafk)

sw$p <- cell_spec(round(sw$p, 4), color =  ifelse(sw$p < 0.05, "red", "white"))

sw %>% kbl(caption = "Outliers identification", escape = F) %>% 
  kable_material_dark("striped")
#No group outliers 

#Teabags 
swt <- Tea %>%
  group_by(Temperature, Treatment, Invertebrates) %>%
  shapiro_test(Teak)

swt$p <- cell_spec(round(swt$p, 4), color =  ifelse(swt$p < 0.05, "red", "white"))

swt %>% kbl(caption = "Outliers identification", escape = F) %>% 
  kable_material_dark("striped")
#no group outliers apart from Inverts yes, t5, control


ggqqplot(Leaf, "Leafk", ggtheme = theme_bw()) +
  facet_grid(Temperature + Treatment ~ Invertebrates, labeller = "label_both")


ggqqplot(Tea, "Teak", ggtheme = theme_bw()) +
  facet_grid(Temperature + Treatment ~ Invertebrates, labeller = "label_both")

#The group outliers shown on this graph ^

##Levenes Test
Leaf %>% levene_test(Leafk ~ Temperature*Treatment*Invertebrates) %>% 
  kbl(caption = "Levene's test") %>% 
  kable_material_dark()

Tea %>% levene_test(Teak ~ Temperature*Treatment*Invertebrates) %>% 
  kbl(caption = "Levene's test") %>% 
  kable_material_dark()
#Levenes Not significant, so we can assume homogeneity of variances in the different groups ^


par(mfrow=c(1,1))
plot(lm(Leafk~Temperature*Treatment*Invertebrates, data=Leaf), 1)
plot(lm(Teak~Temperature*Treatment*Invertebrates,data=Tea),1)
#Some points violate the homogeneity of variances, however the homogeneity has been proved in the Levenes

#Anova
res.aov <- Leaf %>% anova_test(Leafk ~ Temperature*Treatment*Invertebrates)
res.aov %>% kbl(caption = "ANOVA table") %>% 
  kable_material_dark() %>%
  row_spec(c(7), color ="red", bold = T)

res.aov <- Tea %>% anova_test(Teak ~ Temperature*Treatment*Invertebrates)
res.aov %>% kbl(caption = "ANOVA table") %>% 
  kable_material_dark() %>%
  row_spec(c(7), color ="red", bold = T)

resultsleaf <- aov(Leafk ~ Temperature * Treatment * Invertebrates, data = Leaf)
anova(resultsleaf) %>% kbl %>% kable_material_dark()

resultstea <- aov(Teak ~ Temperature * Treatment * Invertebrates, data = Tea)
anova(resultsleaf) %>% kbl %>% kable_material_dark()

report(resultsleaf)
report(resultstea)


#Post Hoc testing 
posthocleaf  <- lm(Leafk ~ Temperature*Treatment*Invertebrates, data = Leaf)
Leaf %>%
  group_by(Temperature) %>%
  anova_test(Leafk ~ Treatment*Invertebrates, error = posthocleaf) %>% 
  kbl(caption = "Two-way interactions") %>% 
  kable_material_dark()

posthoctea  <- lm(Teak ~ Temperature*Treatment*Invertebrates, data = Tea)
Tea %>%
  group_by(Temperature) %>%
  anova_test(Teak ~ Treatment*Invertebrates, error = posthoctea) %>% 
  kbl(caption = "Two-way interactions") %>% 
  kable_material_dark()


#Main Effects 
treatment.effect <- Leaf %>%
  group_by(Temperature, Invertebrates) %>%
  anova_test(Leafk ~ Treatment, error = posthocleaf)

treatment.effect %>% kbl(caption = "Treatment effect") %>% 
  kable_material_dark() %>%
  row_spec(c(1,2), color ="red", bold = T)


treatment.effecttea <- Tea %>%
  group_by(Temperature, Invertebrates) %>%
  anova_test(Teak ~ Treatment, error = posthoctea)

treatment.effecttea %>% kbl(caption = "Treatment effect") %>% 
  kable_material_dark() %>%
  row_spec(c(1,2), color ="red", bold = T)


#Pairwise Comparisons 

pwc15 <- Leaf %>%
  group_by(Temperature, Invertebrates) %>%
  emmeans_test(Leafk ~ Treatment, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk
pwc15 %>% filter(Temperature == "15", Invertebrates == "Yes") %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()

pwleaf <- Leaf %>%
  group_by(Temperature, Invertebrates) %>%
  emmeans_test(Leafk ~ Treatment, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for male at high risk


pwc15 %>% filter(Temperature == "5", Invertebrates == "Yes") %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()


pwleaf
write_csv(pwleaf, "/Users/hamishyoungson/Documents/Imperial/Thesis/Thesis files/Data/Leaf Litter/pwcleaf.csv")


pwctea <- Tea %>%
  group_by(Temperature, Invertebrates) %>%
  emmeans_test(Teak ~ Treatment, p.adjust.method = "bonferroni") %>%
  select(-df, -statistic, -p) # Remove details
# Show comparison results for no invertebrates at 5ºC
pwctea %>% filter(Temperature == "5", Invertebrates == "No ") %>%
  kbl(caption = "Pairwise comparisons") %>% 
  kable_material_dark()

write_csv(pwctea, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Leaf Litter/pwctea.csv")


pwctea





#################### Mortality ###################################


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



Mortsum<- summarySE(Mort,measurevar= "Dead.Gammarus", groupvars=c("Temperature..CºC.", "Treatment"))

Mortsum$Treatment <- factor(Mortsum$Treatment)
ggplot(Mortsum, aes(x= Treatment, y= Dead.Gammarus, fill = factor(Temperature..CºC., levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid", "Fipronil", "P Mix "))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1", "steelblue","lightskyblue")) + 
  labs(fill = "Temperature ºC", title = "") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= Dead.Gammarus-se, ymax= Dead.Gammarus+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Chemical", y="Gammarus Mortality")




bothmort<- aov(Dead.Gammarus ~ Temperature..CºC.*Treatment, data = Mort)
anova(bothmort)
summary(bothmort)

hist(residuals(bothmort))
##mort vs both ^

