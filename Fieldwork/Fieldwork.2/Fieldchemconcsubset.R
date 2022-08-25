rm(list=ls())


library(emmeans)
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
library(rstatix)
library(report)
library(ggpubr)
library(scatterplot3d)


#Read in data
F <- read.csv('Fieldchemconcsubset.csv')
Te <- read.csv('teachemconcsubset.csv')
C <- read.csv('Chlorophyllsubsub.csv')



#Rename
names(F)[34] <- "Urban"
names(Te)[9] <- "Urban"
names(C)[34] <- "Urban"




#Coarse Leaf Litter k calculation 
F$coarseleafk1 <- (log(F$Post.Leaf.Litter.Coarse/F$Pre.Leaf.Litter.Weight..Coarse..g.)/21)*-1

#Fine Leaf Litter k calculation 
F$fineleafk1 <- (log(F$Post.Leaf.Litter.Fine/F$Pre.Leaf.Litter.Weight..Fine..g.)/21)*-1

#teabag 1 k calculation 
Te$teabagk <- (log(Te$Post.Teabag.Weight/Te$Pre.Teabag.Weight..g.)/21)*-1

#wettex fine k calculation (note the log +1!!!)
F$wettexfk <- (log(F$Post.Wettex..Fine/F$Pre.Wettex.Weight)/21)*-1

#wettex coarse k calculation (note the log +1!!!)
F$wettexck <- (log(F$Post.Wettex..Coarse/F$Pre.Wettex.Weight)/21)*-1

#######################Degree Days#############################################

#Degree Day Coarse Leaf Litter k calculation 
F$coarseleafddk1 <- (log(F$Post.Leaf.Litter.Coarse/F$Pre.Leaf.Litter.Weight..Coarse..g.)/(F$Temperature*21))*-1

#Degree Day Fine Leaf Litter k calculation 
F$fineleafddk1 <- (log(F$Post.Leaf.Litter.Fine/F$Pre.Leaf.Litter.Weight..Fine..g.)/(F$Temperature*21))*-1

#Degree Day teabag 1 k calculation 
Te$teabagddk <- (log(Te$Post.Teabag.Weight/Te$Pre.Teabag.Weight..g.)/(Te$Temperature*21))*-1

#Degree Day wettex fine k calculation (note the log +1!!!)
F$wettexfddk <- (log(F$Post.Wettex..Fine/F$Pre.Wettex.Weight)/(F$Temperature*21))*-1

#Degree Day wettex coarse k calculation (note the log +1!!!)
F$wettexcddk <- (log(F$Post.Wettex..Coarse/F$Pre.Wettex.Weight)/(F$Temperature*21))*-1


######degreedays for chlorophyll
C$Chlorophylldd <-  C$Chlorophyll/(C$Temperature*21)

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




#Explanatory variables as factors 
F$Imidacloprid <- factor(F$Imidacloprid)
F$Site <- factor(F$Site)
F$Point <- factor(F$Point)
F$Location <- factor(F$Location)


Te$Imidacloprid <- factor(Te$Imidacloprid)
Te$Site <- factor(Te$Site)
Te$Point <- factor(Te$Point)
Te$Location <- factor(Te$Location)


C$Imidacloprid <- factor(C$Imidacloprid)
C$Site <- factor(C$Site)
C$Point <- factor(C$Point)
C$Location <- factor(C$Location)


########### Linear Mixed Models ###############

################################ LMM's #########################################

#####Coarse Leaf ##############################################################

hist(F$coarseleafddk1)

lmerconcI <- lmer(coarseleafddk1 ~ Imidacloprid.Concentration +(1|Point) , data = F)
lmerconcF <- lmer(coarseleafddk1 ~ Fipronil.Concentration +(1|Point) , data = F)
lmerconcFI <- lmer(coarseleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Point) , data = F)
lmerconcFIT <- lmer(coarseleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Point) , data = F)
lmerconcFIU <- lmer(coarseleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Point) , data = F)

summary(lmerconcI)
summary(lmerconcF)
summary(lmerconcFI)
summary(lmerconcFIT)
summary(lmerconcFIU)


Anova(lmerconcI, test.statistic = "F")
Anova(lmerconcF, test.statistic = "F")
Anova(lmerconcFI, test.statistic = "F")
Anova(lmerconcFIT, test.statistic = "F")
Anova(lmerconcFIU, test.statistic = "F")


# check model assumptions
plot(lmerconcI)
plot(lmerconcF)
plot(lmerconcFI)
plot(lmerconcFIT)
plot(lmerconcFIU)


report(lmerconcI)
report(lmerconcF)
report(lmerconcFI)
report(lmerconcFIT)
report(lmerconcFIU)



#Model of best fit 
model.set <- list(lmerconcI, lmerconcF, lmerconcFI, lmerconcFIT, lmerconcFIU)
model.names <- c("lmerconcI", "lmerconcF" , "lmerconcFI", "lmerconcFIT", "lmerconcFIU")
aictab(model.set, modnames = model.names)




#####Fine Leaf #################################################################


lmerconfI <- lmer(fineleafddk1 ~ Imidacloprid.Concentration +(1|Point) , data = F)
lmerconfF <- lmer(fineleafddk1 ~ Fipronil.Concentration +(1|Point) , data = F)
lmerconfFI <- lmer(fineleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Point) , data = F)
lmerconfFIT <- lmer(fineleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Point) , data = F)
lmerconfFIU <- lmer(fineleafddk1 ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Point) , data = F)

summary(lmerconfI)
summary(lmerconfF)
summary(lmerconfFI)
summary(lmerconfFIT)
summary(lmerconfFIU)


Anova(lmerconfI, test.statistic = "F")
Anova(lmerconfF, test.statistic = "F")
Anova(lmerconfFI, test.statistic = "F")
Anova(lmerconfFIT, test.statistic = "F")
Anova(lmerconfFIU, test.statistic = "F")


# check model assumptions
plot(lmerconfI)
plot(lmerconfF)
plot(lmerconfFI)
plot(lmerconfFIT)
plot(lmerconfFIU)


report(lmerconfI)
report(lmerconfF)
report(lmerconfFI)
report(lmerconfFIT)
report(lmerconfFIU)



#Model of best fit 
model.set <- list(lmerconfI, lmerconfF, lmerconfFI, lmerconfFIT, lmerconfFIU)
model.names <- c("lmerconfI", "lmerconfF" , "lmerconfFI", "lmerconfFIT", "lmerconfFIU")
aictab(model.set, modnames = model.names)




#####Coarse Wettex #################################################################


hist(F$wettexcddk)

F$wettexcddklog <- log(F$wettexcddk) #Logged to increase distributon of normality
hist(F$wettexcddklog)


lmerconwI <- lmer(wettexcddklog ~ Imidacloprid.Concentration +(1|Point) , data = F)
lmerconwF <- lmer(wettexcddklog ~ Fipronil.Concentration +(1|Point) , data = F)
lmerconwFI <- lmer(wettexcddklog ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Point) , data = F)
lmerconwFIT <- lmer(wettexcddklog ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Point) , data = F)
lmerconwFIU <- lmer(wettexcddklog ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Point) , data = F)

summary(lmerconwI)
summary(lmerconwF)
summary(lmerconwFI)
summary(lmerconwFIT)
summary(lmerconwFIU)


Anova(lmerconwI, test.statistic = "F")
Anova(lmerconwF, test.statistic = "F")
Anova(lmerconwFI, test.statistic = "F")
Anova(lmerconwFIT, test.statistic = "F")
Anova(lmerconwFIU, test.statistic = "F")


# check model assumptions
plot(lmerconwI)
plot(lmerconwF)
plot(lmerconwFI)
plot(lmerconwFIT)
plot(lmerconwFIU)


report(lmerconwI)
report(lmerconwF)
report(lmerconwFI)
report(lmerconwFIT)
report(lmerconwFIU)



#Model of best fit 
model.set <- list(lmerconwI, lmerconwF, lmerconwFI, lmerconwFIT, lmerconwFIU)
model.names <- c("lmerconwI", "lmerconwF" , "lmerconwFI", "lmerconwFIT", "lmerconwFIU")
aictab(model.set, modnames = model.names)



##### Fine Wettex ##############################################################


hist(F$wettexfddk)

F$wettexfddklog <- log(F$wettexfddk) #Logged to increase distributon of normality
hist(F$wettexfddklog)

lmerconwfI <- lmer(wettexfddklog ~ Imidacloprid.Concentration +(1|Point) , data = F)
lmerconwfF <- lmer(wettexfddklog ~ Fipronil.Concentration +(1|Point) , data = F)
lmerconwfFI <- lmer(wettexfddklog ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Point) , data = F)
lmerconwfFIT <- lmer(wettexfddklog ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Point) , data = F)
lmerconwfFIU <- lmer(wettexfddklog ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Point) , data = F)

summary(lmerconwfI)
summary(lmerconwfF)
summary(lmerconwfFI)
summary(lmerconwfFIT)
summary(lmerconwfFIU)


Anova(lmerconwfI, test.statistic = "F")
Anova(lmerconwfF, test.statistic = "F")
Anova(lmerconwfFI, test.statistic = "F")
Anova(lmerconwfFIT, test.statistic = "F")
Anova(lmerconwfFIU, test.statistic = "F")


# check model assumptions
plot(lmerconwfI)
plot(lmerconwfF)
plot(lmerconwfFI)
plot(lmerconwfFIT)
plot(lmerconwfFIU)


report(lmerconwfI)
report(lmerconwfF)
report(lmerconwfFI)
report(lmerconwfFIT)
report(lmerconwfFIU)



#Model of best fit 
model.set <- list(lmerconwfI, lmerconwfF, lmerconwfFI, lmerconwfFIT, lmerconwfFIU)
model.names <- c("lmerconwfI", "lmerconwfF" , "lmerconwfFI", "lmerconwfFIT", "lmerconwfFIU")
aictab(model.set, modnames = model.names)


##### Tea ###############################################################


hist(Te$teabagddk)

tlmerconwfI <- lmer(teabagddk ~ Imidacloprid.Concentration +(1|Point) , data = Te)
tlmerconwfF <- lmer(teabagddk~ Fipronil.Concentration +(1|Point) , data = Te)
tlmerconwfFI <- lmer(teabagddk ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Point) , data = Te)
tlmerconwfFIT <- lmer(teabagddk ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Point) , data = Te)
tlmerconwfFIU <- lmer(teabagddk ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Point) , data = Te)

summary(tlmerconwfI)
summary(tlmerconwfF)
summary(tlmerconwfFI)
summary(tlmerconwfFIT)
summary(tlmerconwfFIU)


Anova(tlmerconwfI, test.statistic = "F")
Anova(tlmerconwfF, test.statistic = "F")
Anova(tlmerconwfFI, test.statistic = "F")
Anova(tlmerconwfFIT, test.statistic = "F")
Anova(tlmerconwfFIU, test.statistic = "F")


# check model assumptions
plot(tlmerconwfI)
plot(tlmerconwfF)
plot(tlmerconwfFI)
plot(tlmerconwfFIT)
plot(tlmerconwfFIU)


report(tlmerconwfI)
report(tlmerconwfF)
report(tlmerconwfFI)
report(tlmerconwfFIT)
report(tlmerconwfFIU)



#Model of best fit 
model.set <- list(tlmerconwfI, tlmerconwfF, tlmerconwfFI, tlmerconwfFIT, tlmerconwfFIU)
model.names <- c("tlmerconwfI", "tlmerconwfF" , "tlmerconwfFI", "tlmerconwfFIT", "tlmerconwfFIU")
aictab(model.set, modnames = model.names)


#####Chlorophyll ##############################################################

hist(C$Chlorophylldd)
C$Chlorolog <- log(C$Chlorophylldd) ###Logged to increase normality
hist(C$Chlorolog)

clmerconcI <- lmer(Chlorolog ~ Imidacloprid.Concentration +(1|Location) , data = C)
clmerconcF <- lmer(Chlorolog~ Fipronil.Concentration +(1|Location) , data = C)
clmerconcFI <- lmer(Chlorolog ~ Fipronil.Concentration*Imidacloprid.Concentration +(1|Location) , data = C)
clmerconcFIT <- lmer(Chlorolog ~ Fipronil.Concentration*Imidacloprid.Concentration*Temperature +(1|Location) , data = C)
clmerconcFIU <- lmer(Chlorolog ~ Fipronil.Concentration*Imidacloprid.Concentration*Urban +(1|Location) , data = C)

summary(clmerconcI)
summary(clmerconcF)
summary(clmerconcFI)
summary(clmerconcFIT)
summary(clmerconcFIU)


Anova(clmerconcI, test.statistic = "F")
Anova(clmerconcF, test.statistic = "F")
Anova(clmerconcFI, test.statistic = "F")
Anova(clmerconcFIT, test.statistic = "F")
Anova(clmerconcFIU, test.statistic = "F")

# check model assumptions
plot(clmerconcI)
plot(clmerconcF)
plot(clmerconcFI)
plot(clmerconcFIT)
plot(clmerconcFIU)


report(clmerconcI)
report(clmerconcF)
report(clmerconcFI)
report(clmerconcFIT)
report(clmerconcFIU)



#Model of best fit 
model.set <- list(clmerconcI, clmerconcF, clmerconcFI, clmerconcFIT, clmerconcFIU)
model.names <- c("clmerconcI", "clmerconcF" , "clmerconcFI", "clmerconcFIT", "clmerconcFIU")
aictab(model.set, modnames = model.names)

####PLOTS####################################################################

###Coarse Leaf #############

ggplot(F, aes(x=Imidacloprid.Concentration, y=coarseleafddk1))+
  geom_point() + theme_bw() +
  labs(title = "Coarse Leaf",
       x = "Imidacloprid Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))


ggplot(F, aes(x=Fipronil.Concentration, y=coarseleafddk1))+
  geom_point()  + theme_bw() +
  labs(title = "Coarse Leaf",
       x = "Fipronil Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))






####Fine Leaf    



ggplot(F, aes(x=Imidacloprid.Concentration, y=fineleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "Fine Leaf",
       x = "Imidacloprid Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))


ggplot(F, aes(x=Fipronil.Concentration, y=fineleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Fine Leaf",
       x = "Fipronil Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))

####Coarse Wettex 



ggplot(F, aes(x=Imidacloprid.Concentration, y=wettexcddklog))+
  geom_point() +  theme_bw()  +
  labs(title = "",
       x = "Imidacloprid Concentration",
       y = "Decay Constant k (DD-1)")


ggplot(F, aes(x=Fipronil.Concentration, y=wettexcddklog))+
  geom_point() +  theme_bw() +
  labs(title = "",
       x = "Fipronil Concentration",
       y = "Decay Constant k (DD-1)")


####Fine Wettex 



ggplot(F, aes(x=Imidacloprid.Concentration, y=wettexfddklog))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "",
       x = "Imidacloprid Concentration",
       y = "Decay Constant k (DD-1)")


ggplot(F, aes(x=Fipronil.Concentration, y=wettexfddklog))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Fipronil Concentration",
       y = "Decay Constant k (DD-1)")


####Tea 



ggplot(Te, aes(x=Imidacloprid.Concentration, y=teabagddk))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "Tea",
       x = "Imidacloprid Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))


ggplot(Te, aes(x=Fipronil.Concentration, y=teabagddk))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Tea",
       x = "Fipronil Concentration",
       y= expression(paste("Decay constant k " (DD^-1))))




####Chlorophyll

ggplot(C, aes(x=Imidacloprid.Concentration, y=Chlorolog))+
  geom_point()  + theme_bw()  +
  labs(title = "",
       x = "Imidacloprid Concentration",
       y = "Chlorophyll Production (DD-1)")


ggplot(C, aes(x=Fipronil.Concentration, y=Chlorolog))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Fipronil Concentration",
       y = "Chlorophyll Production (DD-1)")






#########Invertebrates ##################

#Gammarus Imida

F$loggam <- log(F$Gammarus+1)

ggplot(F, aes(x= Imidacloprid.Concentration, y=loggam))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "",
       x = "Imidacloprid Concentration (µg/L)",
       y = "Gammarus abundance (log x+1)")



lmergam <- lmer(loggam ~ Imidacloprid.Concentration +(1|Point) , data = F)

Anova(lmergam, test.statistic = "F")
hist(residuals(lmergam)) ##Gammarus logged to increase normality

report(lmergam)
##Gammarus Fipronil


ggplot(F, aes(x= Fipronil.Concentration, y=loggam))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Fipronil Concentration (µg/L)",
       y = "Gammarus Abundance log (x+1)")


lmergamf <- lmer(Gammarus~ Fipronil.Concentration +(1|Point) , data = F)

Anova(lmergamf, test.statistic = "F")
hist(residuals(lmergamf))


report(lmergamf)

#Asellus
ggplot(F, aes(x=Imidacloprid.Concentration, y=Asellus))+
  geom_point() + theme_bw()  +
  labs(title = "",
       x = "Imidacloprid Concentration",
       y = "Asellus Abundance")


lmeras <- lmer(Asellus ~ Imidacloprid.Concentration +(1|Point) , data = F)

Anova(lmeras, test.statistic = "F")
hist(residuals(lmeras))

report(lmeras)


##Fipronil
ggplot(F, aes(x=Fipronil.Concentration, y=Asellus))+
  geom_point()  + theme_bw() +
  labs(title = "",
       x = "Fipronil Concentration",
       y = "Asellus Abundance")

lmerasf <- lmer(Asellus ~ Fipronil.Concentration +(1|Point) , data = F)

Anova(lmerasf, test.statistic = "F")
hist(residuals(lmerasf))

report(lmerasf)



########Concentration vs Urban 
###Imida
ggplot(F, aes(x=Urban, y=Imidacloprid.Concentration))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "",
       x = "Urban/Suburban landcover (%)",
       y = "Imidacloprid Concentration (µg/L)")


urbanaov <- aov(Imidacloprid.Concentration ~ Urban, data = F)
anova(urbanaov)
plot(urbanaov)
hist(residuals(urbanaov)) ##Assumptions not met 


#######
##Fip

ggplot(F, aes(x=Urban, y= Fipronil.Concentration))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw()  +
  labs(title = "",
       x = "Urban/Suburban landcover (%)",
       y = "Fipronil Concentration (µg/L)")


urbanaovf <- aov(Fipronil.Concentration ~ Urban, data = F)
anova(urbanaovf)
plot(urbanaovf)
hist(residuals(urbanaovf)) ##Assumptions not met 




