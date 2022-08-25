rm(list=ls()) 


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
library(dplyr) 
library(SciViews)
library(ggplot2)
library(dplyr)
library(viridis)
library(wesanderson)
library(gridExtra)
library(tidyverse)
library(ggsci)
library(mvnormtest)
library(IDPmisc)
library(broom)
library(ggpubr)
library(rstatix)
library(report)


# Function gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
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




#read in file
C <- read.csv('chemsubset.csv') 
T <- read.csv("teasubset.csv")
Ch <- read_csv('Chlorophyllsubsub.csv')


C$Imidacloprid <- factor(C$Imidacloprid)

str(C)

#name was far too long
names(C)[33] <- "Urban" 


### k calculation ####

#Coarse Leaf Litter k calculation 
C$coarseleafk1 <- (log(C$Post.Leaf.Litter.Coarse/C$Pre.Leaf.Litter.Weight..Coarse..g.)/21)*-1

#Fine Leaf Litter k calculation 
C$fineleafk1 <- (log(C$Post.Leaf.Litter.Fine/C$Pre.Leaf.Litter.Weight..Fine..g.)/21)*-1

#teabag 1 k calculation 
T$teabagk <- (log(T$Post.Teabag.Weight/Te$Pre.Teabag.Weight..g.)/21)*-1

#wettex fine k calculation (note the log +1!!!)
C$wettexfk <- (log(C$Post.Wettex..Fine/C$Pre.Wettex.Weight)/21)*-1

#wettex coarse k calculation (note the log +1!!!)
C$wettexck <- (log(C$Post.Wettex..Coarse/C$Pre.Wettex.Weight)/21)*-1

#######################Degree Days#############################################

#Degree Day Coarse Leaf Litter k calculation 
C$coarseleafddk1 <- (log(C$Post.Leaf.Litter.Coarse/C$Pre.Leaf.Litter.Weight..Coarse..g.)/(C$Temperature*21))*-1

#Degree Day Fine Leaf Litter k calculation 
C$fineleafddk1 <- (log(C$Post.Leaf.Litter.Fine/C$Pre.Leaf.Litter.Weight..Fine..g.)/(C$Temperature*21))*-1

#Degree Day teabag 1 k calculation 
T$teabagddk <- (log(T$Post.Teabag.Weight/T$Pre.Teabag.Weight..g.)/(T$Temperature*21))*-1

#Degree Day wettex fine k calculation (note the log +1!!!)
C$wettexfddk <- (log(C$Post.Wettex..Fine/C$Pre.Wettex.Weight)/(C$Temperature*21))*-1

#Degree Day wettex coarse k calculation (note the log +1!!!)
C$wettexcddk <- (log(C$Post.Wettex..Coarse/C$Pre.Wettex.Weight)/(C$Temperature*21))*-1


ggboxplot(C, x = "Imidacloprid", y = "coarseleafddk1", col= "blue", ylab= "kDD", main="Coarse Leaf")
ggboxplot(C, x = "Imidacloprid", y = "fineleafrootk", col= "blue",ylab= "kDD" , main="Fine Leaf")
ggboxplot(C, x = "Imidacloprid", y = "finewettexrootk", col="blue", ylab= "kDD", main="Fine Wettex")
ggboxplot(C, x = "Imidacloprid", y = "coarsewettexrootk", col="blue", ylab= "kDD",main="Coarse Wettex" )
ggboxplot(T, x = "Imidacloprid", y = "teabagddk", col= "blue", ylab= "kDD", main="Tea")

C$teabagddk <- (T$teabagddk)


#Diagnostics

#Coarse Leaf
modellc  <- lm(coarseleafddk1 ~ Imidacloprid, data = C)
par(mfrow=c(2,3))
plot(modellc)
hist(modellc$residuals)
shapiro_test(residuals(modellc))

#fine Leaf 
modellf  <- lm((fineleafddk1^(1/3)) ~ Imidacloprid, data = C) #cuberoot to increase normality
par(mfrow=c(2,3))
plot(modellf)
hist(modellf$residuals)
shapiro_test(residuals(modellf))

#Fine Wettex
modelwf  <- lm((wettexfddk^(1/3))  ~ Imidacloprid, data = C) ##Rooted by 1/4 to increase normality 
par(mfrow=c(2,3))
plot(modelwf)
hist(modelwf$residuals)
shapiro_test(residuals(modelwf)) 


#Coarse Wettex 
modelwc  <- lm((wettexcddk^(1/4))  ~ Imidacloprid, data = C)## rooted by 1/4 to increase normality 
par(mfrow=c(2,3))
plot(modelwc)
hist(modelwc$residuals)
shapiro_test(residuals(modelwc))

#Tea
modelt  <- lm(teabagddk  ~ Imidacloprid, data = C)
par(mfrow=c(2,3))
plot(modelt)
hist(modelt$residuals)
shapiro_test(residuals(modelt))


#Anovas
?anova_test

#Anova Test Coarse Leaf
res.aovcl <- C %>% anova_test(coarseleafddk1 ~ Imidacloprid)
res.aovcl


write_csv(res.aovcl, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/rescl.csv")

#Anova Coarse Leaf
resultscl <- aov(coarseleafddk1~ Imidacloprid, data = C)
resultscl
resultscl<- anova(resultscl) 
resultscl

report(resultscl)

#Anova Test Fine Leaf
res.aovfl <- C %>% anova_test((fineleafddk1^(1/3)) ~ Imidacloprid)
res.aovfl

write_csv(res.aovfl,  "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/resfl.csv")

#Anova Fine Leaf
resultsfl <- aov((fineleafddk1^(1/3)) ~ Imidacloprid, data = C)
resultsfl
resultsfl<- anova(resultsfl) 
resultsfl

report(resultsfl)


#Anova Test Fine Wettex
res.aovfw <- C %>% anova_test((wettexfddk^(1/3)) ~ Imidacloprid)
res.aovfw

write_csv(res.aovfw, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/resfw.csv")

#Anova Fine Wettes
resultsfw <- aov((wettexfk^(1/3)) ~ Imidacloprid, data = C)
resultsfw
resultsfw<- anova(resultsfw) 
resultsfw


report(resultsfw)

#Anova Test Coarse Wettex
res.aovcw <- C %>% anova_test((wettexcddk^(1/4)) ~ Imidacloprid)
res.aovcw

write_csv(res.aovcw, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/rescw.csv")
#Anova Coarse Wettex
resultscw <- aov((wettexcddk^(1/4)) ~ Imidacloprid, data = C)
resultscw
resultscw<- anova(resultscw) 
resultscw


report(resultscw)



#Anova Test Tea
res.aovt <- C %>% anova_test(teabagddk~ Imidacloprid)
res.aovt

write_csv(res.aovt, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/rest.csv")

#Anova  Tes
resultst <- aov(teabagddk ~ Imidacloprid, data = C)
resultst
resultst<- anova(resultst) 
resultscw


report(resultscw)


###Barcharts.... 

#Transformation to increase normality 
C$finewettexrootk <- C$wettexfk^(1/3)
C$fineleafrootk <- C$fineleafddk1^(1/3)
C$coarsewettexrootk <- C$wettexcddk^(1/4)


coarseleafsum<- summarySE(C,measurevar= "coarseleafddk1", groupvars=c("Imidacloprid"),na.rm = TRUE)
fineleafsum<- summarySE(C,measurevar= "fineleafrootk", groupvars=c("Imidacloprid"),na.rm = TRUE)
finewettexsum<- summarySE(C,measurevar= "finewettexrootk", groupvars=c("Imidacloprid"),na.rm = TRUE)
coarsewettexsum <- summarySE(C,measurevar= "coarsewettexrootk", groupvars=c("Imidacloprid"),na.rm = TRUE)
teasum <- summarySE(C,measurevar= "teabagddk", groupvars=c("Imidacloprid"),na.rm = TRUE)

names(coarseleafsum)[3] <- "kdd" 
names(fineleafsum)[3] <- "kdd" 
names(finewettexsum)[3] <- "kdd" 
names(coarsewettexsum)[3] <- "kdd" 
names(teasum)[3] <- "kdd" 

coarseleafsum$type <- c(rep("Coarse leaf", 2))
fineleafsum$type <- c(rep("Fine leaf", 2))
finewettexsum$type <- c(rep("Fine Wettex", 2))
coarsewettexsum$type <- c(rep("Coarse Wettex", 2))
teasum$type <- c(rep("Tea", 2))

sum <-as.data.frame(rbind(coarseleafsum,fineleafsum,finewettexsum,coarsewettexsum,teasum)) 


coarseleafsum$Imidacloprid <- factor(coarseleafsum$Imidacloprid)
ggplot(coarseleafsum, aes(x= Imidacloprid, y= kdd, fill= Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Coarse Leaf") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= kdd-se, ymax= kdd+se),  width=.22, position=position_dodge(.9)) +scale_fill_manual(values = c('turquoise','deepskyblue4')) +
  labs (x= "",  y= expression(paste("Decay constant k " (DD^-1))))+ guides(fill = "none")  



fineleafsum$Imidacloprid <- factor(fineleafsum$Imidacloprid)
ggplot(fineleafsum, aes(x= Imidacloprid, y= kdd, fill= Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Fine Leaf") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= kdd-se, ymax= kdd+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "",  y= expression(paste("Decay constant cbrt k " (DD^-1)))) +scale_fill_manual(values = c('turquoise','deepskyblue4'))+ 
  guides(fill = "none")  


finewettexsum$Imidacloprid <- factor(finewettexsum$Imidacloprid)
ggplot(finewettexsum, aes(x= Imidacloprid, y= kdd, fill= Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Fine Wettex") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= kdd-se, ymax= kdd+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "",  y= expression(paste("Decay constant cbrt k " (DD^-1)))) + scale_fill_manual(values = c('turquoise','deepskyblue4')) + 
  guides(fill = "none")  


coarsewettexsum$Imidacloprid <- factor(coarsewettexsum$Imidacloprid)
ggplot(coarsewettexsum, aes(x= Imidacloprid, y= kdd, fill= Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Coarse Wettex") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= kdd-se, ymax= kdd+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "",  y= expression(paste("Decay constant qdrt k " (DD^-1)))) + scale_fill_manual(values = c('turquoise','deepskyblue4'))+ 
  guides(fill = "none")  


teasum$Imidacloprid <- factor(teasum$Imidacloprid)
ggplot(teasum, aes(x= Imidacloprid, y= kdd, fill= Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Tea") +
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= kdd-se, ymax= kdd+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "",  y= expression(paste("Decay constant k " (DD^-1)))) + scale_fill_manual(values = c('turquoise','deepskyblue4')) + 
  guides(fill = "none")  



##Manova

res.man <- manova(cbind(coarseleafddk1, fineleafddk1, wettexfddk,wettexcddk, teabagddk) ~ Imidacloprid, data = C)
summary(res.man)
summary.aov(res.man)


####imidacloprid predicted by land use 

##Binomial GLM Model 
urbanglm <- glm(Imida_bi ~ Urban, family = binomial(link = "logit"), data = C)
uglm <- summary (urbanglm)
uglm <- anova_test(urbanglm)
uglm

report(urbanglm)

write_csv(uglm, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/uglm.csv")

ggboxplot(C, x = "Urban", y = "Imida_bi", col= "blue", ylab= "", main="")
ggboxplot(C, x = "Imida_bi", y = "Urban", col= "black", main="", xlab="Imidacloprid Present", ylab= "Urban/Suburban Landcover (%)")
ggbarplot(C, x = "Imida_bi", y = "Urban", col= "black", main="", xlab="Imidacloprid Present", ylab= "Urban/Suburban Landcover (%)")




confint(urbanglm)
exp(coef(urbanglm))

C$logit <- predict(urbanglm, type = 'link')
C$p <- predict(urbanglm, type = 'response')

modela <- glm(Imida_bi ~ Urban, family = binomial, data = C)
xv <- seq(from = 0, to = 9, by = 0.01)
yv <- predict(object = modela, newdata = list(Percentage.Urban.and.Suburban = xv), type = "response")

par(mfrow = c(1, 1))
plot(C$Urban, C$Imida_bi)
lines(xv, yv)


pred1 <- data.frame(Urban = seq(from = 0,
                                     to = 100, by = 5))
Pred <- predict(urbanglm, newdata = pred1, type = "response")

plot(x = C$Urban, y = C$Imida_bi, ylab= "Imidacloprid", xlab ="Urban/Suburban Landcover (%)")
lines(pred1$Urban, Pred)


#####Linear model k and Urbanisation 


#histo
hist(C$coarseleafddk1)
hist(C$fineleafrootk)
hist(C$finewettexrootk)
hist(C$coarsewettexrootk)
hist(C$teabagddk)
hist(C$Urban)


#plots
plot(coarseleafddk1 ~ Urban, data = C)

plot(fineleafrootk ~ Urban, data = C)

plot(finewettexrootk ~ Urban, data = C)

plot(coarsewettexrootk~Urban, data= C)

plot(teabagddk~Urban, data= C)


#Lm
coarseleaflm<- lm(coarseleafddk1 ~ Urban, data = C)
summary(coarseleaflm)

report(coarseleaflm)
CL <- anova_test(coarseleaflm)

write_csv(CL, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/CL.csv")

###
fineleaflm <- lm(fineleafrootk ~ Urban, data= C)
summary(fineleaflm)

report(fineleaflm)
FL <- anova_test(fineleaflm)
write_csv(FL, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/FL.csv")

###
coarsewettexlm <- lm(coarsewettexrootk ~ Urban, data= C)
summary(coarsewettexlm)

report(coarsewettexlm)
CW <- anova_test(coarsewettexlm)
write_csv(CW, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/CW.csv")

####
finewettexlm <- lm(finewettexrootk ~ Urban, data= C)
summary(finewettexlm)

report(finewettexlm)
FW <- anova_test(finewettexlm)
write_csv(FW, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/FW.csv")

####
tealm <- lm(teabagddk ~ Urban, data= C)
summary(tealm)

report(tealm)
Tea <- anova_test(tealm)
write_csv(Tea, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Subset/Tea.csv")


### Diagnositcs
par(mfrow=c(2,2))
plot(coarseleaflm)
plot(fineleaflm)
plot(coarsewettexlm)
plot(finewettexlm)
plot(tealm)

par(mfrow=c(1,1))



#plots 2 
ggplot(C, aes(x=Urban, y=coarseleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Coarse Leaf",
       x = "Urban and Suburban area (%) ",
       y = "DD")

ggplot(C, aes(x=Urban, y=fineleafrootk))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Fine Leaf",
       x = "Urban and Suburban area (%) ",
       y = "DD")

ggplot(C, aes(x=Urban, y=coarsewettexrootk))+
  geom_point() + theme_bw() +
  labs(title = "Coarse Wettex",
       x = "Urban and Suburban area (%) ",
       y = "DD")

ggplot(C, aes(x=Urban, y=finewettexrootk))+
  geom_point()+ geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Fine Wettex",
       x = "Urban and Suburban area (%) ",
       y = "DD")

ggplot(C, aes(x=Urban, y=teabagddk))+
  geom_point() + theme_bw() +
  labs(title = "Teabag",
       x = "Urban and Suburban area (%) ",
       y = "DD")






