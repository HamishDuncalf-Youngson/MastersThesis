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


#Data 
F <- read.csv("field bioassays2.csv")
Te <- read.csv('teabags.csv')
C <- read.csv('Chlorophyllsub.csv')

#Calculate k values
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

#Rename
names(F)[34] <- "Urban"
names(Te)[9] <- "Urban"
names(C)[34] <- "Urban"



######degreedays for chlorophyll
C$Chlorophylldd <-  C$Chlorophyll/(C$Temperature*21)

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

str(C)

################################ LMM's #########################################

#####Coarse Leaf ##############################################################

hist(F$coarseleafddk1)


lmercl_IMUR <- lmer(coarseleafddk1 ~ Imidacloprid*Urban +(1|Point) , data = F)
lmercl_IMT <- lmer(coarseleafddk1 ~ Imidacloprid*Temperature +(1|Point) , data = F)
lmercl_UT <- lmer(coarseleafddk1 ~ Urban*Temperature +(1|Point) , data = F)
lmercl_IM <- lmer(coarseleafddk1 ~ Imidacloprid +(1|Point) , data = F)
lmercl_T <- lmer(coarseleafddk1 ~ Temperature +(1|Point) , data = F)
lmercl_U <- lmer(coarseleafddk1 ~ Urban +(1|Point) , data = F)
lmercl_UTI <- lmer(coarseleafddk1 ~ Urban*Temperature*Imidacloprid +(1|Point) , data = F)



summary(lmercl_IMUR)
summary(lmercl_IMT)
summary(lmercl_UT)
summary(lmercl_IM)
summary(lmercl_T)
summary(lmercl_U)
summary(lmercl_UTI)



Anova(lmercl_IMUR, test.statistic = "F")
Anova(lmercl_IMT, test.statistic = "F")
Anova(lmercl_UT, test.statistic = "F")
Anova(lmercl_IM, test.statistic = "F")
Anova(lmercl_T, test.statistic = "F")
Anova(lmercl_U, test.statistic = "F")
Anova(lmercl_UTI, test.statistic = "F")

# check model assumptions
plot(lmercl_IMUR)
plot(lmercl_IMT)
plot(lmercl_UT)
plot(lmercl_IM)
plot(lmercl_T)
plot(lmercl_U)
plot(lmercl_UTI)

report(lmercl_IMUR)
report(lmercl_IMT)
report(lmercl_UT)
report(lmercl_IM)
report(lmercl_T)
report(lmercl_U)
report(lmercl_UTI)


#Model of best fit 
model.set <- list(lmercl_IMUR, lmercl_IMT, lmercl_UT, lmercl_IM, lmercl_T, lmercl_U, lmercl_UTI)
model.names <- c("lmercl_IMUR", "lmercl_IMT" , "lmercl_UT", "lmercl_IM", "lmercl_T", "lmercl_U", "lmercl_UTI")
aictab(model.set, modnames = model.names)


lmer_g<- lmer(coarseleafddk1 ~ Imidacloprid + (1|Point), data = F)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")


#####Fine Leaf #################################################################


hist(F$fineleafddk1)



lmerfl_IMUR <- lmer(fineleafddk1 ~ Imidacloprid*Urban +(1|Point) , data = F)
lmerfl_IMT <- lmer(fineleafddk1 ~ Imidacloprid*Temperature +(1|Point) , data = F)
lmerfl_UT <- lmer(fineleafddk1 ~ Urban*Temperature +(1|Point) , data = F)
lmerfl_IM <- lmer(fineleafddk1 ~ Imidacloprid +(1|Point) , data = F)
lmerfl_T <- lmer(fineleafddk1 ~ Temperature +(1|Point) , data = F)
lmerfl_U <- lmer(fineleafddk1 ~ Urban +(1|Point) , data = F)
lmerfl_UTI <- lmer(fineleafddk1 ~ Imidacloprid*Urban*Temperature +(1|Point) , data = F)


summary(lmerfl_IMUR)
summary(lmerfl_IMT)
summary(lmerfl_UT)
summary(lmerfl_IM)
summary(lmerfl_T)
summary(lmerfl_U)
summary(lmerfl_UTI)


Anova(lmerfl_IMUR, test.statistic = "F")
Anova(lmerfl_IMT, test.statistic = "F")
Anova(lmerfl_UT, test.statistic = "F")
Anova(lmerfl_IM, test.statistic = "F")
Anova(lmerfl_T, test.statistic = "F")
Anova(lmerfl_U, test.statistic = "F")
Anova(lmerfl_UTI, test.statistic = "F")


# check model assumptions
plot(lmerfl_IMUR)
plot(lmerfl_IMT)
plot(lmerfl_UT)
plot(lmerfl_IM)
plot(lmerfl_T)
plot(lmerfl_U)
plot(lmerfl_UTI)

report(lmerfl_IMUR)
report(lmerfl_IMT)
report(lmerfl_UT)
report(lmerfl_IM)
report(lmerfl_T)
report(lmerfl_U)
report(lmerfl_UTI)



#Model of best fit 
model.set <- list(lmerfl_IMUR, lmerfl_IMT, lmerfl_UT, lmerfl_IM, lmerfl_T, lmerfl_U, lmerfl_UTI)
model.names <- c("lmerfl_IMUR", "lmerfl_IMT" , "lmerfl_UT", "lmerfl_IM", "lmerfl_T", "lmerfl_U", "lmerfl_UTI")
aictab(model.set, modnames = model.names)


lmer_g<- lmer(fineleafddk1 ~ Imidacloprid + (1|Point), data = F)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")



#####Coarse Wettex #################################################################


hist(F$wettexcddk)

F$wettexcddklog <- log(F$wettexcddk) #Logged to increase distributon of normality
hist(F$wettexcddklog)

lmerwc_IMUR <- lmer(wettexcddklog ~ Imidacloprid*Urban +(1|Point) , data = F)
lmerwc_IMT <- lmer(wettexcddklog ~ Imidacloprid*Temperature +(1|Point) , data = F)
lmerwc_UT <- lmer(wettexcddklog~ Urban*Temperature +(1|Point) , data = F)
lmerwc_IM <- lmer(wettexcddklog ~ Imidacloprid +(1|Point) , data = F)
lmerwc_T <- lmer(wettexcddklog ~ Temperature +(1|Point) , data = F)
lmerwc_U <- lmer(wettexcddklog ~ Urban +(1|Point) , data = F)
lmerwc_UTI<- lmer(wettexcddklog ~ Imidacloprid*Urban*Temperature +(1|Point) , data = F)



summary(lmerwc_IMUR)
summary(lmerwc_IMT)
summary(lmerwc_UT)
summary(lmerwc_IM)
summary(lmerwc_T)
summary(lmerwc_U)
summary(lmerwc_UTI)


Anova(lmerwc_IMUR, test.statistic = "F")
Anova(lmerwc_IMT, test.statistic = "F")
Anova(lmerwc_UT, test.statistic = "F")
Anova(lmerwc_IM, test.statistic = "F")
Anova(lmerwc_T, test.statistic = "F")
Anova(lmerwc_U, test.statistic = "F")
Anova(lmerwc_UTI, test.statistic = "F")



# check model assumptions
plot(lmerwc_IMUR)
plot(lmerwc_IMT)
plot(lmerwc_UT)
plot(lmerwc_IM)
plot(lmerwc_T)
plot(lmerwc_U)
plot(lmerwc_UTI)

report(lmerwc_IMUR)
report(lmerwc_IMT)
report(lmerwc_UT)
report(lmerwc_IM)
report(lmerwc_T)
report(lmerwc_U)
report(lmerwc_UTI)


#Model of best fit 
model.set <- list(lmerwc_IMUR, lmerwc_IMT, lmerwc_UT, lmerwc_IM, lmerwc_T, lmerwc_U, lmerwc_UTI)
model.names <- c("lmerwc_IMUR", "lmerwc_IMT" , "lmerwc_UT", "lmerwc_IM", "lmerwc_T", "lmerwc_U", "lmerwc_UTI")
aictab(model.set, modnames = model.names)


lmer_g<- lmer(wettexcddklog ~ Imidacloprid + (1|Point), data = F)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")




##### Fine Wettex ##############################################################


hist(F$wettexfddk)

F$wettexfddklog <- log(F$wettexfddk) #Logged to increase distributon of normality
hist(F$wettexfddklog)

lmerwf_IMUR <- lmer(wettexfddklog ~ Imidacloprid*Urban +(1|Point) , data = F)
lmerwf_IMT <- lmer(wettexfddklog ~ Imidacloprid*Temperature +(1|Point) , data = F)
lmerwf_UT <- lmer(wettexfddklog~ Urban*Temperature +(1|Point) , data = F)
lmerwf_IM <- lmer(wettexfddklog ~ Imidacloprid +(1|Point) , data = F)
lmerwf_T <- lmer(wettexfddklog ~ Temperature +(1|Point) , data = F)
lmerwf_U <- lmer(wettexfddklog ~ Urban +(1|Point) , data = F)
lmerwf_UTI <- lmer(wettexfddklog~ Urban*Temperature*Imidacloprid +(1|Point) , data = F)

summary(lmerwf_IMUR)
summary(lmerwf_IMT)
summary(lmerwf_UT)
summary(lmerwf_IM)
summary(lmerwf_T)
summary(lmerwf_U)
summary(lmerwf_UTI)

Anova(lmerwf_IMUR, test.statistic = "F")
Anova(lmerwf_IMT, test.statistic = "F")
Anova(lmerwf_UT, test.statistic = "F")
Anova(lmerwf_IM, test.statistic = "F")
Anova(lmerwf_T, test.statistic = "F")
Anova(lmerwf_U, test.statistic = "F")
Anova(lmerwf_UTI, test.statistic = "F")

# check model assumptions
plot(lmerwf_IMUR)
plot(lmerwf_IMT)
plot(lmerwf_UT)
plot(lmerwf_IM)
plot(lmerwf_T)
plot(lmerwf_U)
plot(lmerwf_UTI)

report(lmerwf_IMUR)
report(lmerwf_IMT)
report(lmerwf_UT)
report(lmerwf_IM)
report(lmerwf_T)
report(lmerwf_U)
report(lmerwf_UTI)

#Model of best fit 
model.set <- list(lmerwf_IMUR, lmerwf_IMT, lmerwf_UT, lmerwf_IM, lmerwf_T, lmerwf_U, lmerwf_UTI)
model.names <- c("lmerwf_IMUR", "lmerwf_IMT" , "lmerwf_UT", "lmerwf_IM", "lmerwf_T", "lmerwf_U","lmerwf_UTI" )
aictab(model.set, modnames = model.names)


lmer_g<- lmer(wettexfddklog ~ Imidacloprid + (1|Point), data = F)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")


##### Tea ###############################################################


hist(Te$teabagddk)


lmert_IMUR <- lmer(teabagddk ~ Imidacloprid*Urban +(1|Point) , data = Te)
lmert_IMT <- lmer(teabagddk ~ Imidacloprid*Temperature +(1|Point) , data = Te)
lmert_UT <- lmer(teabagddk~ Urban*Temperature +(1|Point) , data = Te)
lmert_IM <- lmer(teabagddk ~ Imidacloprid +(1|Point) , data = Te)
lmert_T <- lmer(teabagddk ~ Temperature +(1|Point) , data = Te)
lmert_U <- lmer(teabagddk ~ Urban +(1|Point) , data = Te)
lmert_UTI <- lmer(teabagddk~ Urban*Temperature*Imidacloprid +(1|Point) , data = Te)

summary(lmert_IMUR)
summary(lmert_IMT)
summary(lmert_UT)
summary(lmert_IM)
summary(lmert_T)
summary(lmert_U)
summary(lmert_UTI)


Anova(lmert_IMUR, test.statistic = "F")
Anova(lmert_IMT, test.statistic = "F")
Anova(lmert_UT, test.statistic = "F")
Anova(lmert_IM, test.statistic = "F")
Anova(lmert_T, test.statistic = "F")
Anova(lmert_U, test.statistic = "F")
Anova(lmert_UTI, test.statistic = "F")

# check model assumptions
plot(lmert_IMUR)
plot(lmert_IMT)
plot(lmert_UT)
plot(lmert_IM)
plot(lmert_T)
plot(lmert_U)
plot(lmert_UTI)

report(lmert_IMUR)
report(lmert_IMT)
report(lmert_UT)
report(lmert_IM)
report(lmert_T)
report(lmert_U)
report(lmert_UTI)



#Model of best fit 
model.set <- list(lmert_IMUR, lmert_IMT, lmert_UT, lmert_IM, lmert_T, lmert_U, lmert_UTI)
model.names <- c("lmert_IMUR", "lmert_IMT" , "lmert_UT", "lmert_IM", "lmert_T", "lmert_U", "lmert_UTI")
aictab(model.set, modnames = model.names)


lmer_g<- lmer(teabagddk ~ Imidacloprid + (1|Point), data = Te)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")


##### Chlorophyll ###############################################################


hist(C$Chlorophylldd)
C$Chlorolog <- log(C$Chlorophylldd) ###Logged to increase normality
hist(C$Chlorolog)

clmert_IMUR <- lmer(Chlorolog ~ Imidacloprid*Urban +(1|Location) , data = C)
clmert_IMT <- lmer(Chlorolog ~ Imidacloprid*Temperature +(1|Location) , data = C)
clmert_UT <- lmer(Chlorolog~ Urban*Temperature +(1|Location) , data = C)
clmert_IM <- lmer(Chlorolog ~ Imidacloprid +(1|Location) , data = C)
clmert_T <- lmer(Chlorolog ~ Temperature +(1|Location) , data = C)
clmert_U <- lmer(Chlorolog ~ Urban +(1|Location) , data = C)
clmert_UTI <- lmer(Chlorolog~ Urban*Temperature*Imidacloprid +(1|Location) , data = C)

summary(clmert_IMUR)
summary(clmert_IMT)
summary(clmert_UT)
summary(clmert_IM)
summary(clmert_T)
summary(clmert_U)
summary(clmert_UTI)


Anova(clmert_IMUR, test.statistic = "F")
Anova(clmert_IMT, test.statistic = "F")
Anova(clmert_UT, test.statistic = "F")
Anova(clmert_IM, test.statistic = "F")
Anova(clmert_T, test.statistic = "F")
Anova(clmert_U, test.statistic = "F")
Anova(clmert_UTI, test.statistic = "F")

# check model assumptions
plot(clmert_IMUR)
plot(clmert_IMT)
plot(clmert_UT)
plot(clmert_IM)
plot(clmert_T)
plot(clmert_U)


report(clmert_IMUR)
report(clmert_IMT)
report(clmert_UT)
report(clmert_IM)
report(clmert_T)
report(clmert_U)





#Model of best fit 
model.set <- list(clmert_IMUR, clmert_IMT, clmert_UT, clmert_IM, clmert_T, clmert_U, clmert_UTI)
model.names <- c("clmert_IMUR", "clmert_IMT" , "clmert_UT", "clmert_IM", "clmert_T", "clmert_U", "clmert_UTI")
aictab(model.set, modnames = model.names)


lmer_g<- lmer(teabagddk ~ Imidacloprid + (1|Point), data = Te)
em <- emmeans(lmer_g, c("Imidacloprid"))
contrast(em, method = "pairwise")


#### Anova's k ~ Location ################################################ 

ggboxplot(F, x = "Location", y = "coarseleafddk1", color= "Location", palette = c("red", "black", "Blue"))
ggbarplot(F, x = "Location", y = "coarseleafddk1", color= "Location", palette = c("red", "black", "Blue"))
####I've already made these bars, just make some stats 

#Coarse Leaf 
modelclla  <- aov(coarseleafddk1~ Location, data = F)

par(mfrow=c(2,3))
plot(modelclla)
hist(residuals(modelclla))



modelclla <- Anova(modelclla)
modelclla

write_csv(modelclla, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modellcll.csv")

#Fine Leaf
modelflla  <- aov(fineleafddk1~ Location, data = F)

par(mfrow=c(2,3))
plot(modelflla)
hist(residuals(modelflla))


TukeyHSD(modelflla)
modelflla <- Anova(modelflla)
modelflla

write_csv(modelflla, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modellfll.csv")

#Coarse Wettex
modelcwla  <- aov(wettexcddklog~ Location, data = F) 

par(mfrow=c(2,3))
plot(modelcwla) #### Use log as wasn't normal 
hist(residuals(modelcwla))


TukeyHSD(modelcwla)
modelcwla <- Anova(modelcwla)
modelcwla

write_csv(modelcwla, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modellcwl.csv")



#Fine Wettex
modelfwla  <- aov(wettexfddklog~ Location, data = F) 

par(mfrow=c(2,3))
plot(modelfwla) #### Use log as wasn't normal 
hist(residuals(modelfwla))

TukeyHSD(modelfwla)
modelfwla <- Anova(modelfwla)
modelfwla

write_csv(modelfwla, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modellfwl.csv")


#Tea 
modeltea  <- aov(teabagddk~ Location, data = Te) 

par(mfrow=c(2,3))
plot(modeltea) 
hist(residuals(modeltea))

TukeyHSD(modeltea)

modeltea <- Anova(modeltea)
modeltea

write_csv(modeltea, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modeltea.csv")

######Chlorophyll

#Tea 
modelchloro  <- aov(Chlorolog~ Location, data = C) 

par(mfrow=c(2,3))
plot(modelchloro) 
hist(residuals(modelchloro))

TukeyHSD(modelchloro)
modelchloro <- Anova(modelchloro)
modelchloro

write_csv(modelchloro, "/Users/hamishyoungson/Documents/Imperial/Thesis/Data/Fieldwork/Fieldwork.2/Stats/modelchloro.csv")


#######Chlorophyll plot by location 
C$Chlorologr <- C$Chlorolog*-1
  
chlorophyllsum<- summarySE(C,measurevar= "Chlorologr", groupvars=c("Location"),na.rm = TRUE)

chlorophyllsum$Location <- factor(chlorophyllsum$Location)
ggplot(chlorophyllsum, aes(x= Location, y= Chlorologr, fill=Location))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition (Temperature Corrected)") + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= Chlorologr-se, ymax= Chlorologr+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Location",  y="Chlorophyll Production (DD-1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))




#### Anova Temperature ~ Urban Area ##########################################

F$logurban <- log(F$Urban)
F$sqrturban <- sqrt(F$Urban)

ggplot(F, aes(x=logurban, y=Temperature))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban area (log) ",
       y = "Temperature")


modelurbantemp <- aov(Temperature~ logurban, data = F) 

par(mfrow=c(2,3))
plot(modelurbantemp) #### Use log as wasn't normal 
hist(residuals(modelurbantemp))

Anova(modelurbantemp)
report(modelurbantemp)



####### Urban ~ Location ------------

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


urbanlocationsum<- summarySE(F,measurevar= "Urban", groupvars=c("Location"),na.rm = TRUE)

urbanlocationsum$Location <- factor(urbanlocationsum$Location)
ggplot(urbanlocationsum, aes(x= Location, y= Urban, fill=Location)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= Urban-se, ymax= Urban +se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Percentage of urban/suburban landcover (%)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))



####### K ~ Gammarus Abundance (for coarse leaaf and wettex) ####################

#Coarse Leaf
F$loggam <- log(F$Gammarus+1) #logged to increase normality 

lmer_gammarus<- lmer(coarseleafddk1~ loggam+ (1|Point), data = F, na.action = na.omit)
summary(lmer_gammarus)
Anova(lmer_gammarus,  test.statistic = "F")

# check model assumptions
par(mfrow=c(1,1))
plot(lmer_gammarus)

#plot
ggplot(F, aes(x=loggam, y=coarseleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Gammarus Abundance (Log (x+1)) ",
       y= expression(paste("Decay constant k " (DD^-1))))

report(lmer_gammarus)


###########################################

#Coarse Wettex

lmer_gammaruscw<- lmer(wettexcddk~ loggam+ (1|Point), data = F, na.action = na.omit)
summary(lmer_gammaruscw)
Anova(lmer_gammaruscw,  test.statistic = "F")

# check model assumptions
plot(lmer_gammaruscw)

#plot
ggplot(F, aes(x=loggam, y=wettexcddk))+
  geom_point()  + theme_bw() +
  labs(title = "",
       x = "Gammarus Abundance (log +1) ",
       y = "Coarse Wettex kDD")

report(lmer_gammaruscw)



######Impact of Imidacloprid on Gammarus ###########



ImidaGam<- summarySE(F,measurevar= "loggam", groupvars=c("Imidacloprid"),na.rm = TRUE)

ImidaGam$Imidacloprid <- factor(ImidaGam$Imidacloprid)
ggplot(ImidaGam, aes(x= Imidacloprid, y= loggam, fill=Imidacloprid)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= loggam-se, ymax= loggam +se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Imidacloprid",  y="Gammarus Abundance (log x+1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4')) + theme(legend.position = "none") 




lmergambi <- lmer(loggam ~ Imidacloprid +(1|Point) , data = F)
Anova(lmergambi, test.statistic = "F")



plot(lmergambi) #### Use log as wasn't normal 
hist(residuals(lmergambi))


report(lmergambi)


######Impact of Imidacloprid on Asellus ###########

F$logas <- log(F$Asellus+1) ##### Logged to improve normality (slightly?)

ImidaAs<- summarySE(F,measurevar= "logas", groupvars=c("Imidacloprid"),na.rm = TRUE)

ImidaAs$Imidacloprid <- factor(ImidaAs$Imidacloprid)
ggplot(ImidaAs, aes(x= Imidacloprid, y= logas, fill=Imidacloprid)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= logas-se, ymax= logas +se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Imidacloprid",  y="Asellus Abundance (log x+1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4')) + theme(legend.position = "none") 





lmerasbi <- lmer(logas ~ Imidacloprid +(1|Point) , data = F)
Anova(lmerasbi, test.statistic = "F")

 #### Use log as wasn't normal 
hist(residuals(lmerasbi))
plot(lmerasbi)


report(lmerasbi)



#####Leaves ###
leafcoarseddsum<- summarySE(F,measurevar= "coarseleafddk1", groupvars=c("Imidacloprid"),na.rm = TRUE)
leafcoarseddsum2 <- leafcoarseddsum

leaffineddsum<- summarySE(F,measurevar= "fineleafddk1", groupvars=c("Imidacloprid"),na.rm = TRUE)
leaffineddsum2 <- leaffineddsum


names(leafcoarseddsum2)[3] <- "k" 
names(leaffineddsum2)[3] <- "k" 
leafddsum <-as.data.frame(rbind(leafcoarseddsum2,leaffineddsum2)) 
leafddsum$Mesh <- c(rep("coarse", 2),rep("fine",2))

leafddsum$Imidacloprid <- factor(leafddsum$Imidacloprid)
ggplot(leafddsum, aes(x= Imidacloprid, y= k, fill=Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition ") + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= k-se, ymax= k+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Imidacloprid",  y= expression(paste("Decay constant k " (DD^-1)))) + scale_fill_manual(values = c('turquoise','deepskyblue4')) +
  facet_wrap(~Mesh)


######Wettex

#Barchart for Wettex  ####
wettexcoarsesum<- summarySE(F,measurevar= "wettexcddklog", groupvars=c("Imidacloprid"),na.rm = TRUE)
wettexcoarsesum2 <- wettexcoarsesum


wettexfinesum<- summarySE(F,measurevar= "wettexfddklog", groupvars=c("Imidacloprid"),na.rm = TRUE)
wettexfinesum2 <- wettexfinesum


names(wettexfinesum2)[3] <- "k" 
names(wettexcoarsesum2)[3] <- "k" 
wettexsum <-as.data.frame(rbind(wettexcoarsesum2,wettexfinesum2)) 
wettexsum$Mesh <- c(rep("coarse", 2),rep("fine",2))
wettexsum$klog <- wettexsum$k*-1

wettexfinesum <- wettexfinesum
wettexsum$Imidacloprid <- factor(wettexsum$Imidacloprid)
ggplot(wettexsum, aes(x= Imidacloprid, y= klog, fill=Imidacloprid))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Wettex Decomposition") + 
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= klog-se, ymax= klog+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Imidacloprid",   y= expression(paste("Decay constant k " (DD^-1))))   + scale_fill_manual(values = c('turquoise','deepskyblue4'))  +
  facet_wrap(~Mesh)


####Tea 


teabagsum1<- summarySE(Te,measurevar= "teabagddk", groupvars=c("Imidacloprid"),na.rm = TRUE)

teabagsum2 <- teabagsum1
teabagsum2$Imidacloprid <- factor(teabagsum2$Imidacloprid)
ggplot(teabagsum2, aes(x= Imidacloprid, y=teabagddk, fill=Imidacloprid)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Teabag Decomposition") + theme(panel.background = element_rect(fill = NA), 
                                               panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
                                               axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= teabagddk-se, ymax= teabagddk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Imidacloprid",   y= expression(paste("Decay constant k " (DD^-1))))  + scale_fill_manual(values = c('turquoise','deepskyblue4'))



###Chlorophyll
chlorosum1<- summarySE(C,measurevar= "Chlorologr", groupvars=c("Imidacloprid"),na.rm = TRUE)

chlorosum1$Imidacloprid <- factor(chlorosum1$Imidacloprid)
ggplot(chlorosum1, aes(x= Imidacloprid, y=Chlorologr, fill=Imidacloprid)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Chlorophyll") + theme(panel.background = element_rect(fill = NA), 
                                      panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
                                      axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= Chlorologr-se, ymax= Chlorologr+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Imidacloprid",  y= expression(paste("Decay constant k " (DD^-1)))) +  scale_fill_manual(values = c('turquoise','deepskyblue4'))



#####Scatterplot k ~ urban land cover and Temperature
##Coarse Leaf
#plot
ggplot(F, aes(x=Urban, y=coarseleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover",
       y = "Coarse Leaf kDD")


ggplot(F, aes(x=Temperature, y=coarseleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Temperature",
       y = "Coarse Leaf kDD")



F$z_rev <- max(F$Urban) - F$Urban

colors <- c("Blue", "Red", "#56B4E9")
colors <- colors[as.numeric(F$Imidacloprid)]
colors2 <- colors[as.numeric(F$Location)]


s1d<- scatterplot3d(x= F$Temperature, y=F$z_rev, z= F$coarseleafddk1,
              pch =19, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
              main="Coarse Leaf", xlab = "Temperature (ºC)",
              ylab = "Urban/Suburban Landcover (%)",
              zlab = "Decay Constant k (DD-1)")


###Source the function
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

#######Add grids
addgrids3d(x= F$Temperature, y=F$z_rev, z= F$coarseleafddk1, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
        col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)       

                                                                                 

############Fine Leaf 
#plot
ggplot(F, aes(x=Urban, y=fineleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (ºC)",
       y = "Fine Leaf k (DD-1)")


ggplot(F, aes(x=Temperature, y=fineleafddk1))+
  geom_point() +  theme_bw() +
  labs(title = "",
       x = "Temperature (ºC)",
       y = "Fine Leaf k (DD-1)")

### 3D Scatterplot #######


s2d<- scatterplot3d(x= F$Temperature, y=F$z_rev, z= F$fineleafddk1,
                    pch =19, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Fine Leaf", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Decay Constant k (DD-1)")


#######Add grids
addgrids3d(x= F$Temperature, y=F$z_rev, z= F$fineleafddk1, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)       




############ Coarse Wettex ########
#plot
ggplot(F, aes(x=Urban, y=wettexcddklog))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Coarse Wettex k (DD-1)")


ggplot(F, aes(x=Temperature, y=wettexcddklog))+
  geom_point() +  theme_bw() +
  labs(title = "",
       x = "Temperature (ºC)",
       y = "Coarse Wettex k (DD-1)")

### 3D Scatterplot #######
s3d<- scatterplot3d(x= F$Temperature, y=F$z_rev, z= F$wettexcddklog,
                    pch =19, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Coarse Wettex", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Decay Constant k log (DD-1)")


#######Add grids
addgrids3d(x= F$Temperature, y=F$z_rev, z= F$wettexcddklog, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)    



############ Fine Wettex ########
#plot
ggplot(F, aes(x=Urban, y=wettexfddklog))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Fine Wettex k log (DD-1)")


ggplot(F, aes(x=Temperature, y=wettexfddklog))+
  geom_point() +  geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Temperature (ºC)",
       y = "Coarse Wettex k log (DD-1)")

### 3D Scatterplot #######
s3d<- scatterplot3d(x= F$Temperature, y=F$z_rev, z= F$wettexfddklog,
                    pch =19, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Fine Wettex", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Decay Constant k log (DD-1)")


#######Add grids
addgrids3d(x= F$Temperature, y=F$z_rev, z= F$wettexfddklog, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)    



############ Coarse Wettex ########
#plot
ggplot(F, aes(x=Urban, y=wettexcddklog))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Coarse Wettex k (DD-1)")


ggplot(F, aes(x=Temperature, y=wettexcddklog))+
  geom_point() +  theme_bw() +
  labs(title = "",
       x = "Temperature (ºC)",
       y = "Coarse Wettex k (DD-1)")

### 3D Scatterplot #######
s3d<- scatterplot3d(x= F$Temperature, y=F$z_rev, z= F$wettexcddklog,
                    pch =19, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Coarse Wettex", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Decay Constant k log (DD-1)")


#######Add grids
addgrids3d(x= F$Temperature, y=F$z_rev, z= F$wettexcddklog, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)    



############ Tea ########
#plot
ggplot(Te, aes(x=Urban, y=teabagddk))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = " k(DD-1)")


ggplot(Te, aes(x=Temperature, y=teabagddk))+
  geom_point() +  geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Temperature (ºC)",
       y = " k (DD-1)")


Te$z_rev <- max(Te$Urban) - Te$Urban

colors <- c("Blue", "Red")
colors <- colors[as.numeric(Te$Imidacloprid)]



### 3D Scatterplot #######

s3d<- scatterplot3d(x= Te$Temperature, y=Te$z_rev, z= Te$teabagddk,
                    pch =16, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Tea", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Decay Constant k (DD-1)")


#######Add grids
addgrids3d(x= Te$Temperature, y=Te$z_rev, z= Te$teabagddk, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)    


############ Chlorophyll ########
#plot
ggplot(C, aes(x=Urban, y=Chlorologr))+
  geom_point() + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = " Chlorophyll Production log (DD-1)")



ggplot(C, aes(x=Temperature, y=Chlorologr))+
    geom_point() + theme_bw() +
    labs(title = "",
         x = "Temperature (ºC)",
         y = " Chlorophyll Production log (DD-1)")
  

C$z_rev <- max(C$Urban) - C$Urban

colors <- c("Blue", "Red")
colors <- colors[as.numeric(C$Imidacloprid)]



### 3D Scatterplot #######
par(mfrow=c(1,1))
s3d<- scatterplot3d(x= C$Temperature, y=C$z_rev, z= C$Chlorologr,
                    pch =16, color = colors, angle =40, scale.y=1, ylim=c(100, 0),box=FALSE , grid= FALSE, y.ticklabs = c(100, 80, 60, 40, 20, 0), 
                    main="Chlorophyll", xlab = "Temperature (ºC)",
                    ylab = "Urban/Suburban Landcover (%)",
                    zlab = "Chlorophyll Production log (DD-1)")


#######Add grids
addgrids3d(x= C$Temperature, y=C$z_rev, z= C$Chlorologr, grid = c("xy", "xz", "yz"))

####Legend
legend("bottom",s1d$xyz.convert(7.5, 3, 4.5), legend = levels(F$Imidacloprid),
       col =  c("blue", "red"),  pch = 16, inset = -0.32, xpd = TRUE, horiz = TRUE)   


####Gammarus, Asellus ~ Urbanisation #########################


ggplot(F, aes(x=Urban)) + 
  geom_smooth(aes(y = loggam), color = "darkred") + 
  geom_smooth(aes(y = Asellus), color="steelblue", linetype="twodash") 

  
  ggplot(F, aes(x=Urban, y=Gammarus))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Gammarus Abundance")

ggplot(F, aes(x=Urban, y=loggam))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Gammarus Abundance (log)")

ggplot(F, aes(x=Urban, y= Asellus))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Asellus Abundance")

ggplot(F, aes(x=Urban, y= logas))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "",
       x = "Urban/Suburban Landcover (%)",
       y = "Asellus Abundance (log)")



ggplot(F, aes(x=logurban, y = value, color = Invertebrate)) + 
  geom_point(aes(y = loggam , col = "Gammarus (log)")) + 
  geom_smooth(aes(y = loggam , col = "Gammarus (log)"), method="lm") +
  geom_point(aes(y = logas , col = "Asellus (log)")) + theme_bw()+ 
  geom_smooth(aes(y = logas , col = "Asellus (log)"), method="lm") +
labs(title = "",
     x = "Urban/Suburban Landcover (log)",
     y = "Invertebrate Abundance")


##Gammarus Aov
urbanaovgam <- aov(Gammarus ~ Urban, data = F)
anova(urbanaovgam)
plot(urbanaovgam)
hist(residuals(urbanaovgam)) ##Assumptions not met, logged below 

##Gammarus Aov logged
urbanaovgaml <- aov(loggam ~ logurban, data = F)
anova(urbanaovgaml)
plot(urbanaovgaml)
hist(residuals(urbanaovgaml)) 

report(urbanaovgaml)

##Asellus Aov
urbanaovasell <- aov(Asellus ~ Urban, data = F)
anova(urbanaovasell)
plot(urbanaovasell)
hist(residuals(urbanaovasell)) ##Assumptions not met, logged below 

#Asellus Aov logged
urbanaovasella <- aov(logas ~ logurban, data = F)
anova(urbanaovasella)
plot(urbanaovasella)

hist(residuals(urbanaovasella)) 

report(urbanaovasella)



###########################################
###########################################
###########################################
###########################################
#####Change Accordingly...... and add chlorophylll
######### k~ Urban

#####Linear model k and Urbanisation 


#histo
hist(F$coarseleafddk1)
hist(F$fineleafddk1)
hist(F$wettexcddk)
hist(F$wettexfddk)
hist(Te$teabagddk)
hist(C$Chlorophylldd)


#plots
plot(coarseleafddk1 ~ Urban, data = F)

plot(fineleafddk1 ~ Urban, data = F)

plot(wettexfddk ~ Urban, data = F)

plot(wettexcddk~Urban, data= F)

plot(teabagddk~Urban, data= Te)

plot(Chlorophylldd~Urban, data= C)

F$urbanlog <- log(F$Urban)
C$urbanlog <- log(C$Urban)
Te$urbanlog <- log(Te$Urban)



#Lm
coarseleaflm<- lm(coarseleafddk1 ~ Urban, data = F)
summary(coarseleaflm)

report(coarseleaflm)
anova_test(coarseleaflm)



###
fineleaflm <- lm(fineleafddk1 ~ Urban, data= F)
summary(fineleaflm)

report(fineleaflm)
anova_test(fineleaflm)

###
coarsewettexlm <- lm(wettexcddklog ~ Urban, data= F)
summary(coarsewettexlm)

report(coarsewettexlm)
anova_test(coarsewettexlm)


####
finewettexlm <- lm(wettexfddk ~ Urban, data= fine2)
summary(finewettexlm)

report(finewettexlm)
anova_test(finewettexlm)

####
tealm <- lm(teabagddk ~ Urban, data= Te)
summary(tealm)

report(tealm)
anova_test(tealm)


####
chlolm <- lm(Chlorophylldd ~ Urban, data= chloro2)
summary(chlolm)

report(chlolm)
anova_test(chlolm)

### Diagnositcs
par(mfrow=c(2,2))
plot(coarseleaflm)
plot(fineleaflm)
plot(coarsewettexlm) ###logged as wasnt fulfilling assumptions
plot(finewettexlm) #Remove outliers as assumptions not fulfilled. 
plot(tealm)
plot(chlolm) ## Remove Outliers as assumptions not fulfilled. 

chloro2 <- C[-c(1,8), ]
fine2 <- F[-c(5,3), ]

par(mfrow=c(1,1))



#plots 2 
ggplot(F, aes(x=Urban, y=coarseleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Coarse Leaf",
       x = "Urban and Suburban area (%) ",
       y = "k (DD-1)")

ggplot(F, aes(x=Urban, y=fineleafddk1))+
  geom_point() + geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Fine Leaf",
       x = "Urban and Suburban area (%) ",
       y = "k (DD-1)")

ggplot(F, aes(x=Urban, y=(wettexcddklog*-1))) +
  geom_point() + geom_smooth(method="lm", col="black") +theme_bw() +
  labs(title = "Coarse Wettex log",
       x = "Urban and Suburban area (%) ",
       y = "k (DD-1)")

ggplot(fine2, aes(x=Urban, y=wettexfddk))+
  geom_point()+ geom_smooth(method="lm", col="black") + theme_bw() +
  labs(title = "Fine Wettex",
       x = "Urban and Suburban area (%) ",
       y = "k (DD-1)")

ggplot(Te, aes(x=Urban, y=teabagddk))+
  geom_point() + theme_bw() +
  labs(title = "Teabag",
       x = "Urban and Suburban area (%) ",
       y = "k (DD-1)")

ggplot(C, aes(x=Urban, y= Chlorophylldd))+
  geom_point() + theme_bw() +
  labs(title = "Chlorophyll",
       x = "Urban and Suburban area (%) ",
       y = "Chloropyll Production (DD-1)")










