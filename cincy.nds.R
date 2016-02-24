#17-Feb-16
#Script to analyze cincy nds from buried stream study 
#two data files.  cincy.nds is the data file with all the individuals reps for the nds analysis
  #based on prior analysis by Jake, "buried up" and "buried down" are the same, so analyze them as "buried"
  #use cincy.nds to analyze carbon limitation patterns
#cincy.reach contains the averages from the carbon*stream*buried*season nds to study relationships
  #with reach-scale data


############################################################################################
#load and set up data

#packages
install.packages("nlme")
install.packages("nortest")
install.packages("dplyr")
library(nlme)
library(nortest)
library(dplyr)

#attach and evaluate data
nds<-read.table(file="cincy.nds.csv", header=T, sep=',')
str(nds)
unique(nds$stream)
unique(nds$season)
unique(nds$reach)
unique(nds$carbon)

nds$carbon <- factor(nds$carbon, levels = c("Control", "Arabinose", "Cellobiose", "Glucose"))
  #same thing?  nds$carbon <- as.factor(nds$carbon)

############################################################################################
#C limitation analysis with all the data
  #to determine carbon limitation patterns in the aggregate,
    #I think it's better to do it with NRR which takes into account season/reach differences
    #in the response over background.  In that case, the controls are not a value that is 
    #analyzed in the model since the carbon responses are divided by the control.  
    #To continue using cr.area, we can either make a variable or just a find and replace
  #I'll put code and diagnostics for each for now, so we can decide based on results
############################################################################################
fall<-subset(nds, season=="Fall")
spring<-subset(nds, season=="Spring")
fall.spring<-rbind(fall,spring)

M1<-gls(nrr~carbon+season+reach+stream, data=fall.spring, na.action=na.omit)
summary(M1)
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
plot(M1)
hist(residuals(M1))
plot(fall.spring$nrr,residuals(M1)) #not plotting due to NAs in nrr
plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), # code to filter out null value (dplyr package) -J
     residuals(M1), 
     xlab="carbon type", 
     ylab="Residuals")

vf<-varIdent(form=~1|carbon)
  #I tried all effects as an alternate variance structure, and this one produced the best model
M2<-gls(nrr~carbon+season+reach+stream, data=fall.spring, na.action=na.omit, weights=vf)
anova(M1,M2)
summary(M2)
  #arabinose is different than control (p<0.0001), not different than glucose, cellobiose higher (p=0.0426)
  #streams are different, so good blocking factor
  #buried has a stronger response than daylight (p<0.0001)
  #spring a stronger response than fall
qqnorm(residuals(M2))
qqline(residuals(M2))
ad.test(residuals(M2))
plot(M2)
hist(residuals(M2))
plot(fall.spring$nrr,residuals(M2))

vf<-varIdent(form=~1|carbon*season)
M3<-gls(nrr~carbon+season+reach+stream, data=fall.spring, na.action=na.omit, weights=vf)
anova(M1,M2,M3)
  #this has lower AIC, but the conclusion is the same
qqnorm(residuals(M3))
qqline(residuals(M3))
ad.test(residuals(M3))
plot(M3)
hist(residuals(M3))
plot(fall.spring$nrr,residuals(M3))


vf<-varIdent(form=~1|carbon*season*reach)
M4<-gls(nrr~carbon+season+reach+stream, data=fall.spring, na.action=na.omit, weights=vf)
anova(M1,M2,M3,M4)
  #this has even lower AIC, but the conclusion is the same

plot(M4)
hist(residuals(M4))
plot(fall.spring$nrr,residuals(M4))
qqnorm(residuals(M4))

############################################################################################
#C limitation analysis by stream*season*reach
  #this analysis tells us which deployment showed C limitation
############################################################################################

#Amberly
#subset data by stream*season*reach
amberly<-subset(nds,stream=="Amberly")
amberly.fall<-subset(amberly,season=="Fall")
amberly.spring<-subset(amberly,season=="Spring")
amberly.summer<-subset(amberly,season=="Summer")
amberly.fall.buried<-subset(amberly.fall,reach=="buried")
amberly.fall.daylight<-subset(amberly.fall,reach=="daylight")
amberly.spring.buried<-subset(amberly.spring,reach=="buried")
amberly.spring.daylight<-subset(amberly.spring,reach=="daylight")
amberly.summer.buried<-subset(amberly.summer,reach=="buried")
amberly.summer.daylight<-subset(amberly.summer,reach=="daylight")

############################################################################################

#Amberly, fall buried CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.fall.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0001), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.6953

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.fall.buried)
  #variance is unequal among strata p=2.2e-5

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #accounting for strata variance is better

############################################################################################

#Amberly, fall buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.fall.buried, na.action=na.omit)
summary(M1)
  #arabinose is barely lower than the control (p=0.498), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are terrible

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #won't plot because of a null value in cr.afdm
bartlett.test(cr.afdm~carbon, data=amberly.fall.buried)
  #variance is unequal among strata p=3.3e-8

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=amberly.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types
  
anova(M1,M2)
AIC(M1,M2)
  #AIC looks better in the model accounting for strata variance

############################################################################################

#Amberly, fall daylight CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.fall.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is not different than control (p=0.102), not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not great p=0.04343

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.fall.daylight)
  #variance is unequal among strata p=0.00021

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p=0.0091), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model accounting for stratum variance is better

############################################################################################

#Amberly, fall daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.fall.daylight, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.498), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.786

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=amberly.fall.daylight)
  #variance is not different among strata p=0.6403

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=amberly.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)
  
anova(M1,M2)
AIC(M1,M2)
#no significant differences among carbon types
  #AIC looks better in the model without stratum variance and the models are not different

############################################################################################

#Amberly, spring buried CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.spring.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p=0.007), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.3083

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.spring.buried)
  #variance is unequal among strata p=0.0014

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.spring.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#Amberly, spring buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.spring.buried, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.388), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.3135

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error because of two null value
bartlett.test(cr.afdm~carbon, data=amberly.spring.buried)
  #variance is not different among strata p=0.4213

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=amberly.spring.buried, na.action=na.omit,weights=vf)
summary(M2)

anova(M1,M2)
AIC(M1,M2)
  #AIC looks better in the first model, and models are not significantly different

############################################################################################

#Amberly, spring daylight CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.spring.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is not different than control (p=0.11), and not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.2369

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.spring.daylight)
  #variance is not different among strata p=0.0804, but they don't look great

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is not significantly different than control (p=0.0777), and not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #gls model has lower AIC, but models are not significantly different

############################################################################################

#Amberly, spring daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.spring.daylight, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control or glucose (p=0.388), but cellobiose is higher than arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are terrible p=0.001912

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error due to one null value
bartlett.test(cr.afdm~carbon, data=amberly.spring.daylight)
  #variance is different among strata p=0.00072

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=amberly.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance is significantly better 

############################################################################################

#Amberly, summer buried CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.summer.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p=0.004), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.08272

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.summer.buried)
  #variance is unequal among strata p=0.0004

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0002), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#Amberly, summer buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.summer.buried, na.action=na.omit)
summary(M1)
  #arabinose is different than the control (p=0.0451), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are terrible p=4.8e-6

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error because of two null value
bartlett.test(cr.afdm~carbon, data=amberly.summer.buried)
  #variance is different among strata p=0.00085

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=amberly.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences p=0.0697

anova(M1,M2)
AIC(M1,M2)
  #The model with strata variance is significantly better

############################################################################################

#Amberly, summer daylight CR expressed by area

M1<-gls(cr.area~carbon,data=amberly.summer.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p=0.0273), and not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are bad p=0.009145

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.summer.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=amberly.summer.daylight)
  #variance is different among strata p=9.2e-6

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=amberly.summer.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly different than control (p=0.0022) and cellobiose (0.0372),
    #and not different than glucose

anova(M1,M2)
AIC(M1,M2)
#model with strata variance is better

############################################################################################

#Amberly, summer daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=amberly.summer.daylight, na.action=na.omit)
summary(M1)
  #no differences
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
#residuals are terrible p=4.2e-11

hist(residuals(M1), xlab="residuals",main="")
plot(amberly.summer.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error due to three null values
bartlett.test(cr.afdm~carbon, data=amberly.summer.daylight)
  #variance is different among strata p=2.9e-12

M2<-gls(cr.afdm~carbon, data=amberly.summer.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
#model with stratum variance is significantly better 

############################################################################################
#Eastgate
#subset data by stream*season*reach
eastgate<-subset(nds,stream=="Eastgate")
eastgate.fall<-subset(eastgate,season=="Fall")
eastgate.spring<-subset(eastgate,season=="Spring")
eastgate.summer<-subset(eastgate,season=="Summer")
eastgate.fall.buried<-subset(eastgate.fall,reach=="buried")
eastgate.fall.daylight<-subset(eastgate.fall,reach=="daylight")
eastgate.spring.buried<-subset(eastgate.spring,reach=="buried")
eastgate.spring.daylight<-subset(eastgate.spring,reach=="daylight")
eastgate.summer.buried<-subset(eastgate.summer,reach=="buried")
eastgate.summer.daylight<-subset(eastgate.summer,reach=="daylight")

############################################################################################

#eastgate, fall buried CR expressed by area

M1<-gls(cr.area~carbon,data=eastgate.fall.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0001), and also cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.08401

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=eastgate.fall.buried)
  #variance is not different among strata p=0.06866

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=eastgate.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), and also cellobiose and glucose

anova(M1,M2)
AIC(M1,M2)
  #accounting for strata variance is a better AIC, but models not significantly different

############################################################################################

#eastgate, fall buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=eastgate.fall.buried, na.action=na.omit)
summary(M1)
  #no significant differences
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are terrible p=6.7e-5

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #looks a bit weird with glucose
bartlett.test(cr.afdm~carbon, data=eastgate.fall.buried)
  #variance is unequal among strata p=8.3e-5

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=eastgate.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #AIC looks better in the model accounting for strata variance

############################################################################################

#eastgate, fall daylight CR expressed by area

M1<-gls(cr.area~carbon,data=eastgate.fall.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p=0.001), not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not good p=0.00011

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=eastgate.fall.daylight)
  #variance is OK among strata p=0.4182

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=eastgate.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p=0.0045), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
#first model is better on AIC, but not significantly different

############################################################################################

#eastgate, fall daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=eastgate.fall.daylight, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.498), but different than cellobiose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
#residuals are not good p=0.01472

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #null value
bartlett.test(cr.afdm~carbon, data=eastgate.fall.daylight)
  #variance is different among strata p=0.0415

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=eastgate.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #arabinose is not different than the control (p=0.2214), but different than cellobiose

anova(M1,M2)
AIC(M1,M2)
  #stratum variance is a better model
############################################################################################

#eastgate, spring buried CR expressed by area

M1<-gls(cr.area~carbon,data=eastgate.spring.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0001), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are off p=0.04296

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=eastgate.spring.buried)
  #variance is unequal among strata p=0.02216

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=eastgate.spring.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#eastgate, spring buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=eastgate.spring.buried, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.985), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
#residuals are not good p=5.8e-8

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error because of two null values
bartlett.test(cr.afdm~carbon, data=eastgate.spring.buried)
#variance is not different among strata p=0.00060

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=eastgate.spring.buried, na.action=na.omit,weights=vf)
summary(M2)

anova(M1,M2)
AIC(M1,M2)
#AIC better in the stratum variance, models is significantly better

############################################################################################

#eastgate, spring daylight CR expressed by area

M1<-gls(cr.area~carbon,data=eastgate.spring.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p=0.0002), and not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not good p=0.0195

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=eastgate.spring.daylight)
  #variance is different among strata p=0.04897

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=eastgate.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly different than control (p<0.0000), and not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model stratum variance is better

############################################################################################

#eastgate, spring daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=eastgate.spring.daylight, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control or others (p=0.5066)
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
#residuals are OK p=0.141

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=eastgate.spring.daylight)
  #variance is different among strata p=0.04846

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=eastgate.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance is significantly better 

############################################################################################

#eastgate, summer buried CR expressed by area

M1<-gls(cr.area~carbon,data=eastgate.summer.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p=0.0001), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not good p=0.01332

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=eastgate.summer.buried)
  #variance is unequal among strata p=0.0004

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=eastgate.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#eastgate, summer buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=eastgate.summer.buried, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control or others (p=0.1036)
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.2276

hist(residuals(M1), xlab="residuals",main="")
plot(eastgate.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=eastgate.summer.buried)
  #variance is not different among strata p=0.0802

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=eastgate.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #arabinose significantly different than control p=0.0026

anova(M1,M2)
AIC(M1,M2)
  #The model with strata variance has lower AIC but is not significantly better

############################################################################################

#no data for eastgate summer daylight due to high flow



############################################################################################
#Este
#subset data by stream*season*reach
este<-subset(nds,stream=="Este")
este.fall<-subset(este,season=="Fall")
este.spring<-subset(este,season=="Spring")
este.summer<-subset(este,season=="Summer")
este.fall.buried<-subset(este.fall,reach=="buried")
este.fall.daylight<-subset(este.fall,reach=="daylight")
este.spring.buried<-subset(este.spring,reach=="buried")
este.spring.daylight<-subset(este.spring,reach=="daylight")
este.summer.buried<-subset(este.summer,reach=="buried")
este.summer.daylight<-subset(este.summer,reach=="daylight")

############################################################################################

#este, fall buried CR expressed by area

M1<-gls(cr.area~carbon,data=este.fall.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0001), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.1205

hist(residuals(M1), xlab="residuals",main="")
plot(este.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.fall.buried)
  #variance is unequal among strata p=3.6e-7

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #accounting for strata variance is better

############################################################################################

#este, fall buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.fall.buried, na.action=na.omit)
summary(M1)
  #no differences
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not good p=0.02718

hist(residuals(M1), xlab="residuals",main="")
plot(este.fall.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #won't plot because of a null value in cr.afdm
bartlett.test(cr.afdm~carbon, data=este.fall.buried)
  #variance is unequal among strata p=0.048

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=este.fall.buried, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #AIC looks better in the model accounting for strata variance

############################################################################################

#este, fall daylight CR expressed by area

M1<-gls(cr.area~carbon,data=este.fall.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p<0.0001), not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.2452

hist(residuals(M1), xlab="residuals",main="")
plot(este.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.fall.daylight)
  #variance is unequal among strata p=0.000266

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model accounting for stratum variance is better

############################################################################################

#este, fall daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.fall.daylight, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.498), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are not OK p=0.00207

hist(residuals(M1), xlab="residuals",main="")
plot(este.fall.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #no plot due to one null
bartlett.test(cr.afdm~carbon, data=este.fall.daylight)
  #variance is different among strata p=0.0008

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=este.fall.daylight, na.action=na.omit,weights=vf)
summary(M2)

anova(M1,M2)
AIC(M1,M2)
  #strata variance is a better model, but no significant differences
  
############################################################################################

#este, spring buried CR expressed by area

M1<-gls(cr.area~carbon,data=este.spring.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.6384

hist(residuals(M1), xlab="residuals",main="")
plot(este.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.spring.buried)
  #variance is unequal among strata p=0.0048

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.spring.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#este, spring buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.spring.buried, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control (p=0.88), cellobiose and glucose do not differ from arabinose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are acceptable p=0.00062

hist(residuals(M1), xlab="residuals",main="")
plot(este.spring.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
  #plot error because of two null values
bartlett.test(cr.afdm~carbon, data=este.spring.buried)
#variance is not different among strata p=0.098

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=este.spring.buried, na.action=na.omit,weights=vf)
summary(M2)

anova(M1,M2)
AIC(M1,M2)
  #AIC looks better in the second model, and models are not significantly different

############################################################################################

#este, spring daylight CR expressed by area

M1<-gls(cr.area~carbon,data=este.spring.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p=0.0001), and not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
#residuals are acceptable p=0.027

hist(residuals(M1), xlab="residuals",main="")
plot(este.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.spring.daylight)
  #variance is different among strata p=0.00089

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly different than control (p<0.0001), and not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #stratum variance is better

############################################################################################

#este, spring daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.spring.daylight, na.action=na.omit)
summary(M1)
  #arabinose is different than the control (p=0.019), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.7126

hist(residuals(M1), xlab="residuals",main="")
plot(este.spring.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=este.spring.daylight)
  #variance is same among strata p=0.2282

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=este.spring.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #first model has lower AIC, but not difference between them

############################################################################################

#este, summer buried CR expressed by area

M1<-gls(cr.area~carbon,data=este.summer.buried, na.action=na.omit)
summary(M1)
  #Arabinose is significantly greater than control (p<0.0000), but not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are garbage p=1.5e-11

hist(residuals(M1), xlab="residuals",main="")
plot(este.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.summer.buried)
  #variance is unequal among strata p=7.3e-12

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly greater than control (p<0.0001) and glucose, but not different than cellobiose or glucose

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance has lower AIC and is significantly better

############################################################################################

#este, summer buried CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.summer.buried, na.action=na.omit)
summary(M1)
  #arabinose is not different than the control or others (p=0.13)
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are terrible p=0.0007

hist(residuals(M1), xlab="residuals",main="")
plot(este.summer.buried$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=este.summer.buried)
  #variance is different among strata p=0.00078

vf<-varIdent(form=~1|carbon)
M2<-gls(cr.afdm~carbon, data=este.summer.buried, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences p=0.0968

anova(M1,M2)
AIC(M1,M2)
  #The model with strata variance is significantly better

############################################################################################

#este, summer daylight CR expressed by area

M1<-gls(cr.area~carbon,data=este.summer.daylight, na.action=na.omit)
summary(M1)
  #Arabinose is different than control (p<0.0000), and not different than cellobiose or glucose
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.7713

hist(residuals(M1), xlab="residuals",main="")
plot(este.summer.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.area~carbon, data=este.summer.daylight)
  #variance is not different among strata p=0.1057

#analyze data accounting for strata variance
vf<-varIdent(form=~1|carbon)
M2<-gls(cr.area~carbon, data=este.summer.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #Arabinose is significantly different than control (p<0.0000)

anova(M1,M2)
AIC(M1,M2)
  #model with strata variance is better, but not significantly different

############################################################################################

#este, summer daylight CR expressed by afdm

M1<-gls(cr.afdm~carbon,data=este.summer.daylight, na.action=na.omit)
summary(M1)
  #arabinose different than the control p=0.024
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
  #residuals are OK p=0.39

hist(residuals(M1), xlab="residuals",main="")
plot(este.summer.daylight$carbon, residuals(M1), xlab="carbon type", ylab="Residuals")
bartlett.test(cr.afdm~carbon, data=este.summer.daylight)
  #variance is different among strata p=0.009

M2<-gls(cr.afdm~carbon, data=este.summer.daylight, na.action=na.omit,weights=vf)
summary(M2)
  #no significant differences among carbon types

anova(M1,M2)
AIC(M1,M2)
  #model with stratum variance is significantly better 

############################################################################################


#somehow I can use this to plot with na?
#myData[!is.na(myData$Column),]
