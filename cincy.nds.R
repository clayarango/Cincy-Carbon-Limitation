#17-Feb-16
#Hi Clay
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
install.packages("lme4")
library(nlme)
library(nortest)
library(dplyr)
library(lme4)

#attach and evaluate data
nds<-read.table(file="cincy.nds.csv", header=T, sep=',')
str(nds)
unique(nds$stream)
unique(nds$season)
unique(nds$reach)
unique(nds$carbon)

nds$carbon <- factor(nds$carbon, levels = c("Control", "Arabinose", "Cellobiose", "Glucose"))

############################################################################################
#C limitation analysis with all the data

#to determine carbon limitation patterns in the aggregate,
  #I think it's better to do it with NRR which takes into account season/reach differences
  #in the response over background for better comparison.  In that case, the controls are 
  #not a value that is analyzed in the model since the carbon responses are divided by the 
  #control.  
############################################################################################
fall<-subset(nds, season=="Fall")
spring<-subset(nds, season=="Spring")
fall.spring<-rbind(fall,spring)
fall.spring<-droplevels(fall.spring)#removes summer a level bc it has no values


#follow Zuur's method for random effects, p.130-140 
#first model does not include stream
M1<-gls(nrr~carbon*season*reach, data=fall.spring, 
        na.action=na.omit)
anova(M1)
qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))
plot(M1)
hist(residuals(M1))

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     residuals(M1), 
     xlab="Carbon Type", 
     ylab="Residuals")
bartlett.test(nrr~carbon, data=fall.spring)
  #differences among carbon variance (p=0.0029)

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     residuals(M1), 
     xlab="Season", 
     ylab="Residuals")
bartlett.test(nrr~season, data=fall.spring)
  #Fall is obviously much more variable than spring (p<2.2e-16)

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     residuals(M1), 
     xlab="Reach", 
     ylab="Residuals")
bartlett.test(nrr~reach, data=fall.spring)
  #differences between reaches (p<2.2e-16)

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     residuals(M1), 
     xlab="Stream", 
     ylab="Residuals")
bartlett.test(nrr~stream, data=fall.spring)
  #differences among streams (p<2.2e-16)


#My personal preference is to avoid transforming data when possibble (see Zuur
#pages 72, 420-421, 448).  Lets try using the weights and
#random arguments to generate a model that passes the diagnostics.

#use stream as a random effect
M2<-lme(nrr~carbon*season*reach, random=~1|stream, 
        data=fall.spring, na.action=na.omit, method="REML")
anova(M1,M2)

anova(M2)
qqnorm(residuals(M2))
qqline(residuals(M2))
ad.test(residuals(M2)) 
  #better, but no good, p2e-11
plot(M2)  
hist(residuals(M2))

E2<-resid(M2, type="normalized")
F2<-fitted(M2)
op<-par(mfrow=c(2,2), mar=c(4,4,3,2))
plot(x=F2, y=E2, xlab="Fitted Values", ylab="Residuals")

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E2, xlab="Carbon Type", ylab="Residuals")
bartlett.test(nrr~carbon, data=fall.spring)
  #differences among carbon variance (p=0.0029)

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E2, xlab="Season", ylab="Residuals")
bartlett.test(nrr~season, data=fall.spring)
#Fall is still more variable than spring (p<2.2e-16)

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E2, xlab="Reach", ylab="Residuals")
bartlett.test(nrr~reach, data=fall.spring)
  #differences between reaches (p<2.2e-16)
par(op)

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E2, xlab="Stream", ylab="Residuals")
bartlett.test(nrr~stream, data=fall.spring)
  #differences among streams (p<2.2e-16)

#Lets use different variance covariates to improve residual
#distributions
vf1 = varIdent(form = ~ 1|season)  # by season
vf1b = varIdent(form = ~1|season*reach) #alternate varIdent by season*reach interaction
vf2 = varPower(form = ~ fitted(.)) # fitted values as variance covariate (also try varExp, varConstPower)
vf3 = varPower(form = ~ fitted(.)|season)# fitted values as variance covariate (by season?)
vf4 = varExp(form = ~ fitted(.))
vf5 = varExp(form = ~ fitted (.)|season)
vf6 = varConstPower(form = ~ fitted(.))
vf7 = varConstPower(form = ~fitted(.)|season)

M3<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = vf1,
        data=fall.spring, na.action=na.omit, method="REML")

M3b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = vf1b,
        data=fall.spring, na.action=na.omit, method="REML")

M4<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf2),
        data=fall.spring, na.action=na.omit, method="REML")

M4b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf2),
        data=fall.spring, na.action=na.omit, method="REML")

M5<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf3),
        data=fall.spring, na.action=na.omit, method="REML")

M5b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf3),
        data=fall.spring, na.action=na.omit, method="REML")

M6<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf4),
        data=fall.spring, na.action=na.omit, method="REML")

M6b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf4),
        data=fall.spring, na.action=na.omit, method="REML")

M7<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf5),
        data=fall.spring, na.action=na.omit, method="REML")

M7b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf5),
        data=fall.spring, na.action=na.omit, method="REML")

ctrl<-lmeControl(returnObject = TRUE)
  #new object to use in control argument of the next 2 models, force object return when no covergence
    #default maxIter is 50, but setting it to 200 still does not produce convergence
    #opt defaults to "nlmind" but setting default to "optim" does not produce convergence
    #optimMethod defaults to "BFGS" but setting it to "L-BFGS-B" does not produce convergence

M8<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf6),
        data=fall.spring, na.action=na.omit, method="REML",
        control = ctrl)

M8b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf6),
        data=fall.spring, na.action=na.omit, method="REML",
        control = ctrl)

M9<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1, vf7),
        data=fall.spring, na.action=na.omit, method="REML",
        control = ctrl)

M9b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf7),
        data=fall.spring, na.action=na.omit, method="REML",
        control = ctrl)
  
anova(M2, M3, M3b, M4, M4b, M5, M5b, M6, M6b, M7, M7b, M8, M8b, M9, M9b) 
  #model 9 wins, but it didn't converge, so model 7b wins
    #use vf1b and vf5 as alternate variance structures

#graphical analysis of M7b versus M9
E7b<-resid(M7b, type="normalized")
F7b<-fitted(M7b)
E9<-resid(M9, type="normalized")
F9<-fitted(M9)
op<-par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(x=F7b, y=E7b, xlab="Fitted Values", ylab="Residuals", main="Model 7b")
plot(x=F9, y=E9, xlab="Fitted Values", ylab="Residuals", main="Model 9")
  #not well distributed around 0 for either, M7b looks a bit better in the large fitted

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E7b, xlab="Carbon Type", 
     ylab="Residuals", main="Model 7b")

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E9, xlab="Carbon Type", 
     ylab="Residuals", main="Model 9")
  #fewer outliers in M7b

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E7b, xlab="Season", 
     ylab="Residuals", main="Model 7b")

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E9, xlab="Season", 
     ylab="Residuals", main="Model 9")
  #about the same

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E7b, xlab="Reach", 
     ylab="Residuals", main="Model 7b")

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E9, xlab="Reach", 
     ylab="Residuals", main="Model 9")
  #M9 a bit better

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E7b, xlab="Reach", 
     ylab="Residuals", main="Model 7b")

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E9, xlab="Reach", 
     ylab="Residuals", main="Model 9")
  #M9 a bit better, but stream funky in both

par(op)

#Try using using stream as variance covariate
vf8 = varIdent(form = ~1|stream)

M10<-lme(nrr~carbon*season*reach, random=~1|stream, 
         weights = varComb(vf1, vf5, vf8),
         data=fall.spring, na.action=na.omit, method="REML",
         control = ctrl)

M10b<-lme(nrr~carbon*season*reach, random=~1|stream, 
        weights = varComb(vf1b, vf5, vf8), # season*reach, fitted, stream
        data=fall.spring, na.action=na.omit, method="REML")


anova(M2, M3, M3b, M4, M4b, M5, M5b, M6, M6b, M7, M7b, M8, M8b, M9, M9b, M10, M10b)   
  #M10b is the winner so use vf1b, vf5, vf8 as alternate variance structures

#graphical analysis of M10 versus M10b
E10<-resid(M10, type="normalized")
F10<-fitted(M10)
E10b<-resid(M10b, type="normalized")
F10b<-fitted(M10b)

op<-par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(x=F10, y=E10, xlab="Fitted Values", ylab="Residuals", main="Model 10")
plot(x=F10b, y=E10b, xlab="Fitted Values", ylab="Residuals", main="Model 10b")

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E10, xlab="Carbon Type", 
     ylab="Residuals", main="Model 10")

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E10b, xlab="Carbon Type", 
     ylab="Residuals", main="Model 10b")
  #not much difference

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E10, xlab="Season", 
     ylab="Residuals", main="Model 10")

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E10b, xlab="Season", 
     ylab="Residuals", main="Model 10b")
  #M10b might look a bit better

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E10, xlab="Reach", 
     ylab="Residuals", main="Model 10")

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E10b, xlab="Reach", 
     ylab="Residuals", main="Model 10b")
  #M10b better

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E10, xlab="Stream", 
     ylab="Residuals", main="Model 10")

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E10b, xlab="Stream", 
     ylab="Residuals", main="Model 10b")
  #by stream is improved in both, better in M10b

qqnorm(E10, main="Model 10")
qqline(E10)
ad.test(E10) 
  #not normal

qqnorm(E10b, main="Model 10b")
qqline(E10b)
ad.test(E10b) 
  #not normal

h<-hist(E10, main="Model 10")
xfit10<-seq(min(E10),max(E10),length=length(E10))
yfit10<-dnorm(xfit10,mean=mean(E10),sd=sd(E10))
yfit10<-yfit10*diff(h$mids[1:2])*length(E10)
lines(xfit10, yfit10, col="blue", lwd=2)

h<-hist(E10b, main="Model 10b")
xfit10b<-seq(min(E10b),max(E10b),length=length(E10b))
yfit10b<-dnorm(xfit10b,mean=mean(E10b),sd=sd(E10b))
yfit10b<-yfit10b*diff(h$mids[1:2])*length(E10b)
lines(xfit10b, yfit10b, col="blue", lwd=2)

par(op)

anova(M10)
anova(M10b)
  #both yield the same significant factors

# Move on to model selection
M10b.full<-lme(nrr~carbon*season*reach, random=~1|stream, 
              weights = varComb(vf1b, vf5, vf8),
              data=fall.spring, na.action=na.omit, method="ML",
              control = ctrl)

anova(M10b.full)

M10b.a<-update(M10b.full, .~. -carbon:season:reach)
anova(M10b.full, M10b.a)
  #three-way interaction can go because the model is significantly better without it
  #the good model is now M10b.a

M10b.b<-update(M10b.a, .~. -carbon:season)
M10b.c<-update(M10b.a, .~. -carbon:reach)
M10b.d<-update(M10b.a, .~. -season:reach)

anova(M10b.a,M10b.b) #p=0.1537
anova(M10b.a,M10b.c) #p=0.4785
anova(M10b.a,M10b.d) #p=2e-4
  #carbon:reach interaction is the least significant and can be dropped
    #M10b.c is the new full model

M10b.e<-update(M10b.c, .~. -carbon:season)
M10b.f<-update(M10b.c, .~. -season:reach)

anova(M10b.c,M10b.e) #p=0.0355
anova(M10b.c,M10b.f) #p<0.0001
  #carbon:season interaction can probably stay based on p-value, but when we
    #examined other model diagnostics, we decided to take it out since this 
    #p-value is on the edge, and the weight of the evidence suggested it wasn't helping
    #the model much.
  #M10b.e is the new full model

M10b.g<-update(M10b.e, .~. -season:reach)
anova(M10b.g, M10b.e)
  #retain season:reach
    #M10b.e is still the best model

M10b.h<-update(M10b.e, .~. -carbon)
anova(M10.e, M10.h)
  #model is not significantly different without carbon, so it can go

M11.full<-lme(nrr~season*reach,
              random=~1|stream, weights=varComb(vf1b,vf5,vf8),
              method="REML", data=fall.spring, na.action=na.omit,
              control=ctrl)

anova(M11.full)
  #all of these terms are significant

E11<-resid(M11.full, type="normalized")
F11<-fitted(M11.full)
qqnorm(E11)
qqline(E11)
ad.test(E11) 
  #not normal, but not as bad as we started

h<-hist(E11, breaks=16, main="Model 11", xlab="Residuals") #defaulted to too few breaks for my taste
xfit11<-seq(min(E11),max(E11),length=length(E11))
yfit11<-dnorm(xfit11,mean=mean(E11),sd=sd(E11))
yfit11<-yfit11*diff(h$mids[1:2])*length(E11)
lines(xfit11, yfit11, col="blue", lwd=2)

op<-par(mfrow=c(2,2), mar=c(4,4,3,2))
plot(x=F11, y=E11, xlab="Fitted Values", ylab="Residuals")

plot(filter(fall.spring, !is.na(nrr)) %>% select(carbon), 
     E11, xlab="Carbon Type", ylab="Residuals")

plot(filter(fall.spring, !is.na(nrr)) %>% select(season), 
     E11, xlab="Season", ylab="Residuals")

plot(filter(fall.spring, !is.na(nrr)) %>% select(reach), 
     E11, xlab="Reach", ylab="Residuals")
par(op)

plot(filter(fall.spring, !is.na(nrr)) %>% select(stream), 
     E11, xlab="Stream", ylab="Residuals")
  #stream looks pretty good

#interaction plot
op<-par(mfrow=c(1,1), mar=c(4,4,2,1))
x<-fall.spring[complete.cases(fall.spring),]#strips NA rows from data for plotting
x<-droplevels(x)#removes levels with no values, such as summer
with(x, 
     interaction.plot(reach,season,nrr, 
                      ylim=c(0,8),lty=c(1,12),lwd=2,ylab="NRR", 
                      xlab="Reach", trace.label="Season"))
  #fall has stronger NRR than spring, 
  #but in spring, buried have greater response to added carbon
  #whereas in fall, daylight reaches have greater response to added carbon

par(op)

#code here as borrowed from Zuur to look for an additive component, which i
  #don't think is necessary
library(lattice)  
xyplot(E11~nrr|season*reach, data=fall.spring, ylab="Residuals",
       xlab="nrr",
       panel=function(x,y){
         panel.grid(h=-1,v=2)
         panel.points(x,y,col=1)
         panel.loess(x,y,span=0.5,col=1,lwd=2)
       })

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
