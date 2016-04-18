#21-Feb-16
#Adding comment to test git commit
#Script to analyze cincy nds from buried stream study 
#two data files.  cincy.nds is the data file with all the individuals reps for the nds analysis
#based on prior analysis by Jake, "buried up" and "buried down" are the same, so analyze them as "buried"
#use cincy.nds to analyze carbon limitation patterns
#cincy.reach contains the averages from the carbon*stream*buried*season nds to study relationships
#with reach-scale data

#spiralling parameter NAs are replaced with mdl where appropriate


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
reach<-read.table(file="cincy.reach.csv", header=T, sep=',')
unique(reach$stream)
unique(reach$season)
unique(reach$reach)
str(reach)

############################################################################################
#data exploration
#no relationship between glucose.nrr and water chem, transient storage
#no relationship with 2 or 1 station metabolism
#no relationship with N uptake Vf
#no relationship with alg.cm2 or bact.cm2 or peri DM or chla
#no relationship with DHA, A.GALA, B.GALA, B.GLUC, XYLO, NACE, ALA,
  #LEU, CELLO, SULF, PHOS, DOPA, DOPAH2, GLYC, DOPA, CQI, BG, POX,
  #LCI, HIX, BIX, FI, P2H

M1<-gls(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit, method="ML")
anova(M1)
  #CBOM and FBOM are positively related to glucose.nrr, but PERI.DM is not
E1<-(residuals(M1))
qqnorm(E1)
qqline(E1)
ad.test(E1)  
  #residuals are normal
hist(E1, xlab="residuals", main="")
  #not awesome
plot(reach$glucose.nrr, E1, xlab="NRR", ylab="Residuals")
  #seems to have an exponential function

#drop peri.dm from the model since it's not significant and compare the two
M2<-gls(glucose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit, method="ML")
anova(M2)

anova(M1,M2)
  #M2 is better
E2<-(residuals(M2))
qqnorm(E2)
ad.test(E2)
hist(E2,xlab="residuals", main="")
  #not awesome
plot(reach$glucose.nrr, E2, xlab="NRR", ylab="Residuals")
  #about the same

#use stream as a random factor
M3<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream, 
        data=reach, na.action=na.omit, method="ML")
anova(M3)
anova(M2,M3)
  #M3 is marginally better

E3<-(residuals(M3))
qqnorm(E3)
ad.test(E3)
  #normal
hist(E3, xlab="residuals", main="")
  #this has better spread
plot(reach$glucose.nrr, E3, xlab="NRR", ylab="Residuals")
  #still yucky

E3.n<-residuals(M3, type="normalized")
coplot(E3.n~glucose.nrr|stream, data=reach, ylab="Normalized Residuals")
  #stream is heterogenous, but accounted for as a random effect
coplot(E3.n~glucose.nrr|season, data=reach, ylab="Normalized Residuals")
  #season is also heterogenous
coplot(E3.n~glucose.nrr|reach, data=reach, ylab="Normalized Residuals")
  #reach is too

#probably need to account for reach and season somehow

#use alternate variance structures to fix the pattern in the residuals
vf1 = varIdent(form = ~ 1|reach)  
vf2 = varPower(form = ~ fitted(.)) 
vf3 = varPower(form = ~ fitted(.)|reach)
vf4 = varExp(form = ~ fitted(.))
vf5 = varExp(form = ~ fitted (.)|reach)
vf6 = varConstPower(form = ~ fitted(.))
vf7 = varConstPower(form = ~fitted(.)|reach)

ctrl<-lmeControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M4<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit, method="ML")

M5<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M6<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M7<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M8<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M9<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M10<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf7, data=reach, na.action=na.omit,
        control=ctrl, method="ML")
  #produces a singularity using "ML" but not using "REML"

anova(M3,M4)
  #M4 is best, accounts for reach

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M11<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf8, data=reach, na.action=na.omit, method="ML")
  #produces a singularity using "ML" but not using "REML"

M12<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit, method="ML")

M13<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf10, data=reach, na.action=na.omit, method="ML")
  #produces a singularity using "ML" but not using "REML"

anova(M4,M12)
  #M12 is best, accounts for stream

E12<-(residuals(M12))
qqnorm(E12)
ad.test(E12)
  #not normal
hist(E12, xlab="residuals", main="")
  #this looks not great
plot(reach$glucose.nrr, E12, xlab="NRR", ylab="Residuals")
  #still linear


#try a different alternate variance structure with fitted values
vf11 = varPower(form = ~ fitted(.)|glucose.nrr)
vf12 = varPower(form = ~ fitted(.)|CBOM.DM)
vf13 = varPower(form = ~ fitted(.)|FBOM.DM)

M14<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf11, data=reach, na.action=na.omit, method="ML")
M15<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit, method="ML")
M16<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit, method="ML")

anova(M12,M14,M15,M16)
  #M14 is slightly better...why are M14,M15,M16 identical?
  #this is dicey due to the identical fit and the crazy F values


#try combining variance structures
M17<-lme(glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=varComb(vf9,vf11), data=reach, na.action=na.omit, method="ML")
anova(M14,M17)  
  #M14 again
anova(M14)

E4<-(residuals(M4))
qqnorm(E4)
ad.test(E4)
  #not normal
hist(E4, xlab="residuals", main="")
  #big outlier
plot(reach$glucose.nrr, E4, xlab="NRR", ylab="Residuals")
  #these residuals look the best of all, but they aren't great

#in balance, M4 or M12 are the current best, but plot diagnostics for M4 look
  #better than M12

#let's try normalizing the data in case that will pull the outlier down into the cloud
reach$L.glucose.nrr<-log(reach$glucose.nrr+1)

#start basic
M18<-lme(L.glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         data=reach, na.action=na.omit, method="ML")
#try with the same weights as the best model above
M19<-lme(L.glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf1, data=reach,na.action=na.omit, method="ML")

AIC(M4,M18,M19)
    #M18 has lower AIC     
anova(M18)

E18<-(residuals(M18))
qqnorm(E18)
qqline(E18)
ad.test(E18)  
#residuals are normal
hist(E18, xlab="residuals", main="")
  #looks OK
plot(reach$L.glucose.nrr, E18, xlab="NRR", ylab="Residuals")
  #these residuals have the best vertical spread of all

anova(M14,M18,M19)#can't use this bc response variables are different
anova(M18)
  #same significant factors
summary(M18)
  #interpretation is that the biofilm response to glucose is strongest with 
    #greater standing stocks of CBOM and FBOM, possibly because of high CBOM
    #standing stocks in the fall

##############################################################################
##############################################################################


M1<-gls(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit)
anova(M1)
  #CBOM significant, FBOM is weak, PERI is not related
E1<-(residuals(M1))
qqnorm(E1)
qqline(E1)
ad.test(E1)  
  #residuals are normal
hist(E1, xlab="residuals", main="")
  #OKish
plot(reach$arabinose.nrr, E1, xlab="NRR", ylab="Residuals")
  #seems to have an exponential function

#drop peri.dm from the model since it's not significant and compare the two
M2<-gls(arabinose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit)
anova(M2)

anova(M1,M2)
#M2 is better, and FBOM is significant without PERI.DM

E2<-(residuals(M2))
qqnorm(E2)
ad.test(E2)
hist(E2,xlab="residuals", main="")
  #OKish
plot(reach$arabinose.nrr, E2, xlab="NRR", ylab="Residuals")
  #about the same

#use stream as a random factor
M3<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream, 
        data=reach, na.action=na.omit, method="REML")
anova(M3)
anova(M2,M3)
  #M3 is marginally better on AIC

E3<-(residuals(M3))
qqnorm(E3)
ad.test(E3)
  #normal
hist(E3, xlab="residuals", main="")
  #this looks funkier
plot(reach$arabinose.nrr, E3, xlab="NRR", ylab="Residuals")
  #still yucky

E3.n<-residuals(M3, type="normalized")
coplot(E3.n~arabinose.nrr|stream, data=reach, ylab="Normalized Residuals")
  #stream is heterogenous, but accounted for as a random effect
coplot(E3.n~arabinose.nrr|season, data=reach, ylab="Normalized Residuals")
  #season is also heterogenous
coplot(E3.n~arabinose.nrr|reach, data=reach, ylab="Normalized Residuals")
  #reach is too

#probably need to account for reach and season somehow

vf1 = varIdent(form = ~ 1|reach)  
vf2 = varPower(form = ~ fitted(.)) 
vf3 = varPower(form = ~ fitted(.)|reach)
vf4 = varExp(form = ~ fitted(.))
vf5 = varExp(form = ~ fitted (.)|reach)
vf6 = varConstPower(form = ~ fitted(.))
vf7 = varConstPower(form = ~fitted(.)|reach)

ctrl<-lmeControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M4<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit)

M5<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit,
        control=ctrl)

M6<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit,
        control=ctrl)
#won't converge

M7<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit,
        control=ctrl)

M8<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit,
        control=ctrl)
#won't converge

M9<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit,
        control=ctrl)

M10<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf7, data=reach, na.action=na.omit,
         control=ctrl)
#won't converge

anova(M2,M3,M4,M5,M7,M9)
anova(M3,M5)
anova(M5)
summary(M5)
  #M5 is best, but no longer significant

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M11<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf8, data=reach, na.action=na.omit)

M12<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit)

M13<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf10, data=reach, na.action=na.omit)

anova(M5,M11,M12,M13)
  #M5 is still best, M12 is nearly as good but it retains CBOM.DM

E5<-(residuals(M5))
qqnorm(E5)
ad.test(E5)
  #not normal
hist(E5, xlab="residuals", main="")
  #big outlier
plot(reach$arabinose.nrr, E5, xlab="NRR", ylab="Residuals")
  #linear

#try a different alternate variance structure with fitted values
vf11 = varPower(form = ~ fitted(.)|arabinose.nrr)
vf12 = varPower(form = ~ fitted(.)|CBOM.DM)
vf13 = varPower(form = ~ fitted(.)|FBOM.DM)

M14<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf11, data=reach, na.action=na.omit)
M15<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit)
M16<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit)

anova(M5,M14,M15,M16)
  #M5 is still best

#try combining variance structures
M17<-lme(arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=varComb(vf2,vf9), data=reach, na.action=na.omit)
anova(M5,M17)
  #M17 better on AIC, but not significantly different

E17<-(residuals(M17))
qqnorm(E17)
ad.test(E17)
  #not normal
hist(E17, xlab="residuals", main="")
  #yuck
plot(reach$arabinose.nrr, E17, xlab="NRR", ylab="Residuals")
  #residuals are still linear or a power function
anova(M17)


#let's try normalizing the data in case that will pull the outlier down into the cloud
reach$L.arabinose.nrr<-log(reach$arabinose.nrr+1)

#start simple
M18<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         data=reach, na.action=na.omit)
#try same weights from the best model above
M19<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=varComb(vf2,vf9), data=reach,na.action=na.omit)
AIC(M17,M18,M19)
  #M18 has lower AIC     

E18<-(residuals(M18))
qqnorm(E18)
qqline(E18)
ad.test(E18)  
  #residuals are normal
hist(E18, xlab="residuals", main="")
  #looks OK
plot(reach$L.arabinose.nrr, E18, xlab="NRR", ylab="Residuals")
  #is this better than the unnormalized?  still has 4 large points, but the vertical spread looks better

anova(M17,M18,M19)#can't use this bc response variables are different
anova(M18)
  #the factors are significant in the simple model with log normalized data
summary(M18)
  #CBOM.DM is not significant with the summary command?

#interpretation is that the biofilm response to arabinose is higher with greater standing stocks
  #of CBOM or FBOM, but the effect size looks very small

##############################################################################
##############################################################################


M1<-gls(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit)
anova(M1)  
E1<-(residuals(M1))
qqnorm(E1)
qqline(E1)
ad.test(E1)  
  #residuals are normal
hist(E1, xlab="residuals", main="")
  #OKish
plot(reach$cellobiose.nrr, E1, xlab="NRR", ylab="Residuals")
  #same yucky pattern we've been seeing

#drop peri.dm from the model since it's not significant and compare the two
M2<-gls(cellobiose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit)
anova(M2)

anova(M1,M2)
#M2 is better, and FBOM can be dropped
M3<-gls(cellobiose.nrr~CBOM.DM,data=reach,na.action=na.omit)
anova(M3)

anova(M2,M3)
#M3 is better

E3<-(residuals(M3))
qqnorm(E3)
qqline(E3)
ad.test(E3)  
  #residuals are normal
hist(E3, xlab="residuals", main="")
  #OKish
plot(reach$cellobiose.nrr, E3, xlab="NRR", ylab="Residuals")
  #exponential relationship

#use stream as a random factor
M4<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream, 
        data=reach, na.action=na.omit, method="REML")
anova(M4)
anova(M3,M4)
#M4 is marginally better on AIC

E4<-(residuals(M4))
qqnorm(E4)
ad.test(E4)
  #normal
hist(E4, xlab="residuals", main="")
  #this looks good
plot(reach$cellobiose.nrr, E4, xlab="NRR", ylab="Residuals")
  #still yucky

E3.n<-residuals(M3, type="normalized")
coplot(E3.n~arabinose.nrr|stream, data=reach, ylab="Normalized Residuals")
  #stream is heterogenous, but accounted for as a random effect
coplot(E3.n~arabinose.nrr|season, data=reach, ylab="Normalized Residuals")
  #season is also heterogenous
coplot(E3.n~arabinose.nrr|reach, data=reach, ylab="Normalized Residuals")
  #reach is too

#probably need to account for reach and season somehow

#try alternate variance structures
vf1 = varIdent(form = ~ 1|reach)  
vf2 = varPower(form = ~ fitted(.)) 
vf3 = varPower(form = ~ fitted(.)|reach)
vf4 = varExp(form = ~ fitted(.))
vf5 = varExp(form = ~ fitted (.)|reach)
vf6 = varConstPower(form = ~ fitted(.))
vf7 = varConstPower(form = ~fitted(.)|reach)

ctrl<-lmeControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M5<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit)

M6<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit,
        control=ctrl)
#singular

M7<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit,
        control=ctrl)
#error

M8<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit,
        control=ctrl)
#won't converge

M9<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit,
        control=ctrl)
#error

M10<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit,
        control=ctrl)

M11<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf7, data=reach, na.action=na.omit,
         control=ctrl)

anova(M4,M5,M10,M11)
  #M4 is still the best

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M12<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf8, data=reach, na.action=na.omit)

M13<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit)

M14<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf10, data=reach, na.action=na.omit)

anova(M4,M12,M13,M14)
anova(M4,M13)
  #M13 is best, accounts for stream as a random factor and as an Ident

E13<-(residuals(M13))
qqnorm(E13)
ad.test(E13)
  #not normal
hist(E13, xlab="residuals", main="")
  #yuck
plot(reach$cellobiose.nrr, E13, xlab="NRR", ylab="Residuals")
  #linear...why are these so linear?  this seems odd

#try a different alternate variance structure with fitted values
vf11 = varPower(form = ~ fitted(.)|cellobiose.nrr)
vf12 = varPower(form = ~ fitted(.)|CBOM.DM)
vf13 = varPower(form = ~ fitted(.)|FBOM.DM)

M15<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf11, data=reach, na.action=na.omit)
M16<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit)
M17<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit)

anova(M13,M15,M16,M17)
  #M15 is best

#try combining variance structures
M18<-lme(cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=varComb(vf9,vf11), data=reach, na.action=na.omit)
anova(M15,M18)
  #M15 better on AIC, but not significantly different

E15<-(residuals(M15))
qqnorm(E15)
ad.test(E15)
  #not normal
hist(E15, xlab="residuals", main="")
  #yuck
plot(reach$cellobiose.nrr, E15, xlab="NRR", ylab="Residuals")
  #residuals are still linear

#let's try normalizing the data
reach$L.cellobiose.nrr<-log(reach$cellobiose.nrr+1)

#start simple
M19<-lme(L.cellobiose.nrr~CBOM.DM, random=~1|stream,
         data=reach, na.action=na.omit)
#try the weights from the best model above
M20<-lme(L.cellobiose.nrr~CBOM.DM, random=~1|stream,
         weights=vf11, data=reach,na.action=na.omit)
AIC(M15,M19,M20)
#M20 has lower AIC     

E20<-(residuals(M20))
qqnorm(E20)
qqline(E20)
ad.test(E20)  
  #residuals are not normal
hist(E20, xlab="residuals", main="")
  #looks OK
plot(reach$L.cellobiose.nrr, E20, xlab="NRR", ylab="Residuals")
  #residuals are still messed up

anova(M15,M19,M20)#can't use this bc response variables are different
anova(M20)
#CBOM is still significant
summary(M20)

#interpretation is that the biofilm response to cellobiose is stronger with lower 
  #standing stocks of CBOM, but the model residuals are screwed up so maybe not valid?
#try to average all carbon sources together because the NDS data showed no carbon effect
  #or do a MANOVA approach to analyze the carbon source simultaneously

##############################################################################
##############################################################################



x<-gls(glucose.nrr~daily.par,data=reach,na.action=na.omit)
anova(x)  
x<-gls(arabinose.nrr~daily.par,data=reach,na.action=na.omit)
anova(x)  
x<-gls(cellobiose.nrr~daily.par,data=reach,na.action=na.omit)
anova(x)
reach$daily.par
  #weak positive relationship between glucose daily PAR, p<0.10 with arabinose, p<0.15 cellobiose
    #i don't think this is a great rabbit hole to go down.  daily.par is a lot of NA and the buried
    #reaches are all zero.  there are just 4 values that aren't zero or NA
    #we've already incorporated reach into some of the models above, so we
    #know reach matters

####################################################################
####################################################################

#N acquistion enzymes

M1<-gls(NACE.DM~season+reach+stream, data=reach,na.action=na.omit, method="ML")
  #can't do interactions bc of just one point per season*reach*stream combination
anova(M1)
summary(M1)
  #significantly lower in spring compared to Fall

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #residuals are normal
plot(M1)
  #looks OK
hist(residuals(M1))
  #not bad
plot(reach$NACE.DM,residuals(M1))
  #positive relationship?
  y<-lm(residuals(M1)~reach$NACE.DM)
  summary(y)
  #yep
    
plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #OK
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK

#looks like stream is weirdest. let's try it as a random factor
M2<-lme(NACE.DM~season+reach, random=~1|stream, 
        data=reach,na.action=na.omit, method="ML")

anova(M1,M2)
  #M1 is better, so maybe nesting is not the best here

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varIdent(form=~1|stream*season*reach)#this one won't work because there are 0 df for the three interaction

M3<-gls(NACE.DM~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit, method="ML")
M4<-gls(NACE.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit, method="ML")
M5<-gls(NACE.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit, method="ML")
M6<-gls(NACE.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit, method="ML")
  #false convergence
M7<-gls(NACE.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit, method="ML")
  #false convergence

anova(M1,M3,M4,M5)#accounting for stream variance is best, which makes sense
  #M3 is best


qqnorm(residuals(M14))
qqline(residuals(M14))
ad.test(residuals(M14))  
  #residuals are normal
plot(M14)
  #not bad, but there's an odd gap between 400-900
hist(residuals(M14))
plot(reach$NACE.DM,residuals(M14))
  #looks OK for the n

plot(reach$season,residuals(M14), xlab="Season", ylab="Residuals")
#OK
plot(reach$reach,residuals(M14), xlab="Reach", ylab="Residuals")
#OK
plot(reach$stream,residuals(M14), xlab="Stream", ylab="Residuals")
#OK

#by stream is still funky so try stream as a random effect along with the varIdent
M8<-lme(NACE.DM~season+reach+stream, random=~1|stream,
        weights=vf1, data=reach,na.action=na.omit)

anova(M3,M8)
  #M3 still wins

#try a model based on fitted values
vf7 = varPower(form = ~ fitted(.))
vf8 = varExp(form = ~ fitted(.))
vf9 = varConstPower(form = ~ fitted(.))

M9<-gls(NACE.DM~season+reach+stream, 
        weights=varComb(vf1,vf7), data=reach,na.action=na.omit, method="ML")

M10<-gls(NACE.DM~season+reach+stream, 
        weights=varComb(vf1,vf8), data=reach,na.action=na.omit, method="ML")
  #doesn't converge

M11<-gls(NACE.DM~season+reach+stream, 
        weights=varComb(vf1,vf9), data=reach,na.action=na.omit, method="ML")

anova(M3,M9,M11)
  #model 9 is the best, but graphical diagnostics look terrible.  stick with M3 

anova(M3)

#term deletion to optimize model, remove stream which is least significant
M12<-gls(NACE.DM~season+reach, 
        weights=vf1, data=reach, na.action=na.omit, method="ML")

anova(M3,M12)
  #model 12 is better than M3 because no significant difference
anova(M12)

M13<-gls(NACE.DM~season, weights=vf1, 
         data=reach, na.action=na.omit, method="ML")

anova(M12,M13)
  #model 3 is the best
M14<-gls(NACE.DM~season, weights=vf1, 
         data=reach, na.action=na.omit, method="REML")

anova(M14)
summary(M14)
  #interpretation is that Fall has the most N acquiring enzymes NACE compared to spring and summer
    #possibly due to a pulse of recalcitrant terrestrial C?  
  #M14 has worse looking diagnostics, but we deleted the insignificant terms 
    #and the interpretation is the same

####################################################################
####################################################################

#N acquisition

M1<-gls(LEU.C~season+reach+stream, data=reach,na.action=na.omit, method="ML")
anova(M1)
summary(M1)
  #almost higher in spring compared to fall p=0.0503

qqnorm(residuals(M5))
qqline(residuals(M5))
ad.test(residuals(M5))  
  #OK, but they look like they have an exponential pattern
plot(M5)
  #variation proportional to fit
hist(residuals(M5))
plot(reach$LEU.C,residuals(M5))
  #positive relationship?
  y<-lm(residuals(M1)~reach$LEU.C)
  summary(y)
    #by outlier...need to deal with this

#try to use stream as a random factor
M2<-lme(LEU.C~season+reach, random=~1|stream, 
          data=reach,na.action=na.omit, method="ML")
anova(M1,M2)
  #M1 is better, so the random factor is not the best

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(LEU.C~season+reach+stream, 
       weights=vf1, data=reach,na.action=na.omit, method="ML")  
M4<-gls(LEU.C~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit, method="ML")  
M5<-gls(LEU.C~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit, method="ML")
M6<-gls(LEU.C~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit, method="ML")
  #error
M7<-gls(LEU.C~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit, method="ML")
  #error
M8<-gls(LEU.C~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit, method="ML")
  #no covergence
M9<-gls(LEU.C~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        method="ML")
  #no covergence

anova(M1,M3,M4,M5)
  #M5 is best, using season as an alternate variance structure

#try a combined variance structure just for fun
M10<-gls(LEU.C~season+reach+stream, 
        weights=varComb(vf3,vf7), data=reach,na.action=na.omit,
        method="ML")  

anova(M5,M10)
  #M8 wins

anova(M5)
summary(M5)
  #no significant factors...Spring is greater than Fall at p=0.0834
    #other similar models have significant factors...ie., M6
    #we should look at the residual details and see if another model 
    #accounts for residuals better

plot(reach$season,residuals(M8), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M8), reach$season)  
  #not OK
plot(reach$reach,residuals(M8), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M8), reach$reach)  
  #OK
plot(reach$stream,residuals(M8), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M8), reach$stream)  
  #OK

qqnorm(residuals(M8))
qqline(residuals(M8))
ad.test(residuals(M8))  
  #they don't look great
plot(M8)
  #these look better
hist(residuals(M8))
plot(reach$LEU.C,residuals(M8))
  #positive relationship?
  y<-lm(residuals(M8)~reach$LEU.C)
  summary(y)
    #driven by the same outlier

#Stream is significant in the model, but this isn't interesting
#We'll also need to do a find+replace to swap LEU.C with LEU.DM to check
  #that variable.  Looks like LEU.C in that some models are significant
  #but others are not

##################################################################
##################################################################

#Sulfur acquisition  

M1<-gls(SULF.DM~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #higher in spring than fall or summer, same with SULF.C

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #OK, but they don't look great
plot(M1)
  #gack
hist(residuals(M1))
plot(reach$SULF.DM,residuals(M1))
  #three major outliers

plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #OK
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK, but those three big outliers might be from Amberly
    #probably need to account for stream in the variance structure

#try to use stream as a random factor
M2<-lme(SULF.DM~season+reach, random=~1|stream, 
        data=reach,na.action=na.omit, method="REML")
anova(M1,M2)
  #M1 is better, so the random factor is not the best

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(SULF.DM~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(SULF.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(SULF.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(SULF.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(SULF.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(SULF.DM~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(SULF.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M6,M7,M8,M9)
  #M9 is best, using fitted values as an alternate variance structure

#try a combined variance structure just for fun
M10<-gls(SULF.DM~season+reach+stream, 
         weights=varComb(vf3,vf7), data=reach,na.action=na.omit)  

anova(M1,M9,M10)
  #M9 wins

anova(M9)
summary(M9)
  #season not significant with anova, but spring is higher than fall with summary

plot(reach$season,residuals(M9), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M9), reach$season)  
  #not OK
plot(reach$reach,residuals(M9), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M9), reach$reach)  
  #OK, but not pretty
plot(reach$stream,residuals(M9), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M9), reach$stream)  
  #OK, but not pretty

qqnorm(residuals(M9))
qqline(residuals(M9))
ad.test(residuals(M9))  
  #the 3 highest and lowest extreme values are freaking out
plot(M9)
  #the variance structure created a bi-modal distribution here
hist(residuals(M9))
plot(reach$SULF.DM,residuals(M9))
  #cloud of most points with the three highest and lowest along the line

#we'll need to explore the individual models in more detail to find the one 
  #that handles the residuals the best...the AIC is misleading us here, 
  #also find and replace SULF.DM with SULF.C

##################################################################
##################################################################

#phosphorus acquisition

M1<-gls(PHOS.DM~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #lower in summer compared to fall, anova p values aren't consistent

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #OK
plot(M1)
  #OK
hist(residuals(M1))
plot(reach$PHOS.DM,residuals(M1))
  #two outliers or increasing trend?

  y<-lm(residuals(M1)~reach$PHOS.DM)
  summary(y)
    #significant,  driven by outliers

plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #OK
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK, but not the best
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK, outliers in Amberly again?

#try to use stream as a random factor
M2<-lme(PHOS.DM~season+reach, random=~1|stream, 
          data=reach,na.action=na.omit, method="REML")
anova(M1,M2)
  #M1 is better, so the random factor is not the best
  
#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))
  
ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)
  
M3<-gls(PHOS.DM~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(PHOS.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(PHOS.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(PHOS.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(PHOS.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(PHOS.DM~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
  #no convergence
M9<-gls(PHOS.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)
  #no convergence

anova(M1,M3,M4,M5,M6,M7)
#M7 is best, using stream*season as an alternate variance structure
  
#try a combined variance structure just for fun
M10<-gls(PHOS.DM~season+reach+stream, 
         weights=varComb(vf1,vf5), data=reach,na.action=na.omit)  
  
anova(M1,M7,M10)
  #M7 wins
  
anova(M7)
summary(M7)
  #those are pretty sick F values...
  
plot(reach$season,residuals(M7), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M7), reach$season)  
  #OK
plot(reach$reach,residuals(M7), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M7), reach$reach)  
  #OK
plot(reach$stream,residuals(M7), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M7), reach$stream)  
  #not great
  
qqnorm(residuals(M7))
qqline(residuals(M7))
ad.test(residuals(M7))  
  #not good, but it could be worse
plot(M7)
  #this looks OK
hist(residuals(M7))
plot(reach$PHOS.DM,residuals(M7))
  #looks better when we started, but those two outliers are still odd

#we'll need to explore the individual models in more detail to find the one 
  #that handles the residuals the best 
  #also find and replace PHOS.DM with PHOS.C

##################################################################
##################################################################

#this is a peroxidase assay that correlates to lignin degradation, so 
  #it's a metric of recalcitrant carbon use

M1<-gls(DOPAH2.DM~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #lower in daylight compared to buried

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #not so good
plot(M1)
  #sort of a U shape
hist(residuals(M1))
  #major outliers
plot(reach$DOPAH2.DM,residuals(M1))
  #increasing trend?

  y<-lm(residuals(M1)~reach$DOPAH2.DM)
  summary(y)
    #significant

plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #OK, but spring is quite variable
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #Este is all over the map, and big outliers in the others

#try to use stream as a random factor
M2<-lme(DOPAH2.DM~season+reach, random=~1|stream, 
          data=reach,na.action=na.omit, method="REML")
anova(M1,M2)
  #M1 is better, so the random factor is not the best
  
#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))
  
ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
  #false convergence
M8<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
  #no convergence
M9<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M6,M9)
  #M1 is best on AIC
  #graphical analysis indicates M4 or M5 might handle residual variation better 

#try a combined variance structure just for fun
M10<-gls(DOPAH2.DM~season+reach+stream, 
         weights=varComb(vf2,vf7), data=reach,na.action=na.omit)  
anova(M1,M4,M5,M10)

plot(reach$season,residuals(M11), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M11), reach$season)  
  #OK
plot(reach$reach,residuals(M11), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M11), reach$reach)  
  #OK
plot(reach$stream,residuals(M11), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M11), reach$stream)  
  #not great
  
qqnorm(residuals(M11))
qqline(residuals(M11))
ad.test(residuals(M11))  
  #not good, but it could be worse
plot(M11)
  #this looks OK
hist(residuals(M11))
plot(reach$DOPAH2.DM,residuals(M11))

  #M10 does some things better but the residuals still scale linearly with the predictor 
    #deleting terms retains reach as significant, but the models look worse
  #can repeat this on DOPA.DM and DOPA.C (phenol oxidase assay...recalcitrant carbon)
    #and with DOPAH2.C, but none of these start significant (some are 0.05<x<0.10)

##################################################################
##################################################################
#GLYC and CQI (index of carbon quality as ratio of GLYC/DOPA, higher numbers are more labile)
  #neither are significant with first model run

##################################################################
##################################################################
#POX is the phenoloxidase assay, so it's a metric of recalcitrant C
  #Brian Hill's alternate calculation, equivalent to Sinsabaugh's DOPA?

M1<-gls(POX~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #almost lower in daylight compared to buried

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #not so good

plot(M1)
  #not horrible
hist(residuals(M1))
  #major outliers
plot(reach$POX,residuals(M1))
  #outliers


plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #not OK big variation in spring
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK

#season looks like the biggest problem
#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(POX~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(POX~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(POX~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(POX~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(POX~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
  #false convergence
M8<-gls(POX~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(POX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M2,M3,M4,M5,M6,M8,M9)
  #not much difference among the models

#try a combined variance structure just for fun
M10<-gls(POX~season+reach+stream, 
         weights=varComb(vf3,vf7), data=reach, na.action=na.omit) 

anova(M1,M5,M10)
  #the combined model is not an improvement

qqnorm(residuals(M11))
qqline(residuals(M11))
ad.test(residuals(M11))  
  #not so good

plot(M11)
  #not horrible
hist(residuals(M11))
  #major outliers
plot(reach$POX,residuals(M11))
  #outliers

plot(reach$season,residuals(M11), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M11), reach$season)  
  #not OK big variation in spring
plot(reach$reach,residuals(M11), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M11), reach$reach)  
  #OK
plot(reach$stream,residuals(M11), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M11), reach$stream)  
  #OK

#try normalizing?
reach$L.POX<-log(reach$POX)

M11<-gls(L.POX~season+reach+stream, data=reach,na.action=na.omit)
anova(M11)
summary(M11)
  #graphical analysis shows a big improvement


##################################################################
##################################################################
#LCI is an index of carbon lability calculcated as recalcitrant/(labile+recalcitrant)
  #so a smaller number means there is more labile carbon and a bigger number means more
  #recalcitrant carbon

M1<-gls(LCI~season+reach+stream, data=reach, na.action=na.omit, method="ML")
anova(M1)
summary(M1)
  #smaller number in daylight compared to buried (more labile in daylight),
    #but the seasonal pattern is opposite what we might expect.
    #this is a ratio rather than a raw number

qqnorm(residuals(M6))
qqline(residuals(M6))
ad.test(residuals(M6))  
  #OK

plot(M6)
  #unimodal or not bad?
hist(residuals(M6))
plot(reach$LCI,residuals(M6))
  #increasing variation with LCI

plot(reach$season,residuals(M6), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M6), reach$season)  
  #looks great
plot(reach$reach,residuals(M6), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M6), reach$reach)  
  #looks great
plot(reach$stream,residuals(M6), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M6), reach$stream)  
  #OK

#mostly looks good except for how the residuals increase with LCI  
#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(LCI~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(LCI~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(LCI~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(LCI~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(LCI~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
  #false convergence
M8<-gls(LCI~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(LCI~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M6,M8,M9)
  #M6 which is stream*reach is best
anova(M6)
summary(M6)
  #model interpretation is better...daylight and spring have lower values (greater carbon lability)
    #than buried and fall

##################################################################
##################################################################
#HIX is from Pennino's work, it is the humification index, so a bigger number
  #means more humic compounds

M1<-gls(HIX~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #almost lower in spring compared to fall, p=0.0538

qqnorm(residuals(M3))
qqline(residuals(M3))
  #my goodness
ad.test(residuals(M3))  
  #how is this not worse than p=0.03?

plot(M3)
  #one major outlier that is underestimated
hist(residuals(M3))
plot(reach$HIX,residuals(M3))
  #one distinct outlier

plot(reach$season,residuals(M3), xlab="Season", ylab="Residuals")
  #there's that outlier
bartlett.test(residuals(M3), reach$season)  
  #OK
plot(reach$reach,residuals(M3), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M3), reach$reach)  
  #not OK
plot(reach$stream,residuals(M3), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M3), reach$stream)  
  #OK, eastgate, buried spring with the outlier
    #this value is in row 44 of Pennino's summary excel file, and it has other
    #unusual values compared to the other samples and seasons.

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(HIX~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(HIX~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(HIX~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(HIX~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(HIX~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(HIX~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(HIX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M6,M7,M8,M9)
  #M3 wins, stream
anova(M3)
  #graphical diagnostics don't show an improvement with the outlier.  we might have to dump it?

##################################################################
##################################################################
#BIX is from Pennino's work, it is the biological freshness index.  Bigger number
#means more autochthony

M1<-gls(BIX~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #spring and summer have higher BIX than fall

qqnorm(residuals(M5))
qqline(residuals(M5))
ad.test(residuals(M5))  
  #good
hist(residuals(M5))
plot(reach$BIX,residuals(M5))
  #major categorical variation...I think I see where Fall is

plot(reach$season,residuals(M5), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M5), reach$season)  
  #OK
plot(reach$reach,residuals(M5), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M5), reach$reach)  
  #OK
plot(reach$stream,residuals(M5), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M5), reach$stream)  
  #OK

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(BIX~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(BIX~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(BIX~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(BIX~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
  #false convergence
M7<-gls(BIX~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
  #false convergence
M8<-gls(BIX~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(BIX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M8,M9)
  #M8 wins, but graphical analysis shows it doesn't handle season well
  #check M5, which had the season varIdent, but it doesn't do it any better
  #is this as good as it gets given the bimodality?
  
##################################################################
##################################################################
#FI is from Pennino's work, it is the fluorescence index.  Bigger number
#bigger=microbial, lower=terrestrial

M1<-gls(FI~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #higher in spring and summer compared to fall

qqnorm(residuals(M9))
qqline(residuals(M9))
ad.test(residuals(M9))  
hist(residuals(M9))
  #big outlier
plot(reach$FI,residuals(M9))
  #major categorical variation

plot(reach$season,residuals(M9), xlab="Season", ylab="Residuals")
  #egad
bartlett.test(residuals(M9), reach$season)  
  #OK on p value, but ugly
plot(reach$reach,residuals(M9), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M9), reach$reach)  
  #OK
plot(reach$stream,residuals(M9), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M9), reach$stream)  
  #OK

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(FI~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(FI~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(FI~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(FI~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
  #false convergence
M7<-gls(FI~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(FI~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(FI~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M7,M8,M9)
  #looks like M9 wins.

#try a combined variance structure just for fun
M10<-gls(FI~season+reach+stream, 
         weights=varComb(vf3,vf7), data=reach, na.action=na.omit) 
anova(M9,M10)
  #M9 still the best
 
  #random thought...is high fall variance compared to narrow summer variance
    #accounted for by different disturbance frq (more t-storm summer) or
    #differences in OM (leaf input in fall)? Maybe we need a covariate like days 
    #since last storm or CBOM standing stock?

##################################################################
##################################################################
#PRO, FUL, HUM are from Pennino's work but he says it's not the best to use the
  #raw numbers and to use the P2H instead.  I'm not workign on those, but I'm leaving
  #this old code here for now

x<-gls(PRO~season+reach+stream, data=reach,na.action=na.omit)#protein value...ratio to H is better per Pennino
summary(x)
  #higher in spring and summer compared to fall

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #outlier
plot(x)
  #outlier
hist(residuals(x))
  #i feel like i'm repeating myself
plot(reach$PRO,residuals(x))
  #major outlier in the predictor

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK on p value, but not good
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK, but not the best
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #eastgate is whack
  #need to look at spring buried eastgate

vf<-varIdent(form=(~1|season))
x2<-gls(PRO~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2)
  #the model with seasonal variance accounted for is better, but need diagnostics


x<-gls(FUL~season+reach+stream, data=reach,na.action=na.omit)
summary(x)
#higher in spring and summer compared to fall and in eastgate compared to Amberly

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #on the edge, but they don't look good
plot(x)
  #not terrible but a possible "U" shaped pattern
hist(residuals(x))
  #lots of underestimates
plot(reach$FUL,residuals(x))
  #not bad

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK on p value, but spring has some wide values
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #eastgate is whack
  
vf<-varIdent(form=(~1|season))
x2<-gls(FUL~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2)
  #the first model is better



x<-gls(HUM~season+reach+stream, data=reach,na.action=na.omit)
summary(x)  
  #higher in spring and summer compared to fall and in eastgate compared to Amberly

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #yuck
plot(x)
  #bimodal estimates...too high or too little
hist(residuals(x))
plot(reach$HUM,residuals(x))
  #not good

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #lots of underestimates

vf<-varIdent(form=(~1|stream))
x2<-gls(HUM~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2)
  #improvement in AIC, but no significant difference

##################################################################
##################################################################
#P2H is the PRO/HUM ratio, bigger numbers are more protein like OM compared to 
  #humic-like OM

M1<-gls(P2H~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)  
  #higher in spring and summer compared to fall

qqnorm(residuals(M7))
qqline(residuals(M7))
ad.test(residuals(M7))  
  #that one outlier is bad, but otherwise OK
hist(residuals(M7))
plot(reach$P2H,residuals(M7))
  #not awful...just that one outlier

plot(reach$season,residuals(M7), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M7), reach$season)  
  #not good
plot(reach$reach,residuals(M7), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M7), reach$reach)  
  #OK
plot(reach$stream,residuals(M7), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M7), reach$stream)  
  #spring buried eastgate, which was the same bit outlier earlier

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(P2H~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(P2H~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(P2H~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(P2H~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)  
M7<-gls(P2H~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(P2H~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit,
        control=ctrl)
M9<-gls(P2H~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit,
        control=ctrl)

anova(M1,M3,M4,M5,M6,M7,M8,M9)
  #M7 look like it wins, stream*season
anova(M7)
  #I don't trust it when the F ratios are so enormous

##################################################################
##################################################################
#CQI and LCI should be inversely related to each other
  #bigger CQI=better carbon, bigger LCI=worse carbon

M1<-gls(CQI~LCI, data=reach, na.action=na.omit)
anova(M1)
summary(M1)  
  #negative relationship

plot(reach$CQI,reach$LCI)
#gack...need to get rid of that outlier...Eastgate Fall Daylight

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #that one outlier is bad, but otherwise OK
plot(x)
  #very distinct pattern due to that big outlier
hist(residuals(x))
  #not good



