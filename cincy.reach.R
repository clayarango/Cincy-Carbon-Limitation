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

M1<-gls(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit, 
        method="REML")  # Zuur says start with REML
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

#Jake: finish optimizing the random part before optimizing the fixed part.
#Jake: comment out variable selection
#drop peri.dm from the model since it's not significant and compare the two
# M2<-gls(glucose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit, method="ML")
# anova(M2)
# 
# anova(M1,M2)
#   #M2 is better
# E2<-(residuals(M2))
# qqnorm(E2)
# ad.test(E2)
# hist(E2,xlab="residuals", main="")
#   #not awesome
# plot(reach$glucose.nrr, E2, xlab="NRR", ylab="Residuals")
#   #about the same


# Jake: continue working on random effects.
#use stream as a random factor
M2<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream, 
        data=reach, na.action=na.omit, 
        method="REML") # use REML for random structure
anova(M2)
anova(M1,M2)
  #M2 is marginally better

E2<-(residuals(M2))
qqnorm(E2)
ad.test(E2)
  #normal
hist(E2, xlab="residuals", main="")
  #this has better spread
plot(reach$glucose.nrr, E2, xlab="NRR", ylab="Residuals")
  #still yucky

E2.n<-residuals(M2, type="normalized")
coplot(E2.n~glucose.nrr|stream, data=reach, ylab="Normalized Residuals")
  #stream is heterogenous, but accounted for as a random effect
coplot(E2.n~glucose.nrr|season, data=reach, ylab="Normalized Residuals")
  #season is also heterogenous
coplot(E2.n~glucose.nrr|reach, data=reach, ylab="Normalized Residuals")
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

M3<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit)

M4<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M5<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit,
        control=ctrl)
  #won't converge

M6<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit)

M7<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit)

M8<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit)
  #won't converge, even when using control = ctrl.  Omitted here to save time.

M9<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf7, data=reach, na.action=na.omit, control=ctrl)
  #produces Error

anova(M1,M2,M3,M6,M7)
  #M6 is best, accounts residuals as a function of fitted values

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M10<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf8, data=reach, na.action=na.omit, 
        method="REML")  # I think method should be REML

M11<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit, 
         method="REML")

M12<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf10, data=reach, na.action=na.omit, 
         method="REML")

anova(M1,M2,M3,M6,M7,M10,M11,M12)
  #M6 is still the best

E6<-(residuals(M6))
qqnorm(E6)
ad.test(E6)
  #not normal
hist(E6, xlab="residuals", main="")
  #this looks not great
plot(reach$glucose.nrr, E6, xlab="NRR", ylab="Residuals")
  #still linear


#try a different alternate variance structure with fitted values
vf11 = varPower(form = ~ fitted(.)|glucose.nrr)
vf12 = varPower(form = ~ fitted(.)|CBOM.DM)
vf13 = varPower(form = ~ fitted(.)|FBOM.DM)

M13<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf11, data=reach, na.action=na.omit, method="REML")
M14<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit, method="REML")
M15<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit, method="REML")

anova(M1,M2,M3,M6,M7,M10,M11,M12,M13,M14,M15)
#M6 is still the best
  #M14 is slightly better...why are M13,M14,M15 identical?
  #this is dicey due to the identical fit and the crazy F values


#try combining variance structures
M16<-lme(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=varComb(vf9,vf11), data=reach, na.action=na.omit, method="REML")
anova(M1,M2,M3,M6,M7,M10,M11,M12,M13,M14,M15,M16)  
  #M6 again
anova(M6)  # nothing significant!

#Does liklihood ratio test indicate that M6 (random + weights) is better than M1 or M2?
anova(M1, M6) #yes
anova(M2, M6) #yes
# See lines 162-166 for model diagnostics.  Look terrible.
# Dont' think we can resolve with weights and random


#let's try normalizing the data in case that will fix the issues
#Zuur pg. 91, step 6 "consider a transformation on the response
#variable as a last resort"
reach$L.glucose.nrr<-log(reach$glucose.nrr+1)

#start basic
M17<-gls(L.glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM,
         data=reach, na.action=na.omit, method="REML")
# Add random
M18<-lme(L.glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         data=reach, na.action=na.omit, method="REML")
#try with the same weights as the best model above
M19<-lme(L.glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf4, data=reach,na.action=na.omit, method="REML")

AIC(M17,M18,M19)
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


# Now it is time to find the optimal fixed structure (pg.91, step 7)
# F-statistic, bad because it uses sequential testing and the order of the main
# effects is important.  Uses REML
anova(M18) # cant trust results

#t-statistic  This works unless one variable is a factor with more than 2 levels (uses REML)
summary(M18) #  PERI.DM not significant, CBOM and FBOM are.

# Remove PERI
M19 <- lme(L.glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
           data=reach,na.action=na.omit, method="REML")

summary(M19) #FBOM and CBOM significant.  Stop here!

#For kicks, lets try optimizing the fixed components using liklihood ratio test
#This requires ML.
#Compare nested models 
M21 <- lme(L.glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream, # full model
           data=reach, na.action=na.omit, method="ML")

M22 <- lme(L.glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream, # remove PERI
           data=reach, na.action=na.omit, method="ML")

M23 <- lme(L.glucose.nrr~CBOM.DM+PERI.DM, random=~1|stream, # remove FBOM
           data=reach, na.action=na.omit, method="ML")

M24 <- lme(L.glucose.nrr~FBOM.DM+PERI.DM, random=~1|stream, # remove CBOM
           data=reach, na.action=na.omit, method="ML")

anova(M21, M22)  # p=0.44 for PERI
anova(M21, M23)  # p=0.01 for FBOM
anova(M21, M24)  # p=0.006 for CBOM

# Remove PERI, least significant term
# New full model = M22
M25 <- lme(L.glucose.nrr~FBOM.DM, random=~1|stream, # remove CBOM
                  data=reach, na.action=na.omit, method="ML")

M26 <- lme(L.glucose.nrr~CBOM.DM, random=~1|stream, # remove FBOM
           data=reach, na.action=na.omit, method="ML")

anova(M22, M25) # CBOM is significant, p = 0.0065
anova(M22, M26) # FBOM is significant, p = 0.0084
# End of Backwards selection to identify optimal fixed structure

# Fit final model
MFinal <- lme(L.glucose.nrr~CBOM.DM+FBOM.DM, random=~1|stream, # remove PERI
              data=reach, na.action=na.omit, 
              method="REML")  # now return to REML!
summary(MFinal)  # same results as summary(M19), line 246 above!

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

#Commented out to optimize random structure before fixed structure

#drop peri.dm from the model since it's not significant and compare the two
#M2<-gls(arabinose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit)
#anova(M2)
#
#anova(M1,M2)
##M2 is better, and FBOM is significant without PERI.DM
#
#E2<-(residuals(M2))
#qqnorm(E2)
#ad.test(E2)
#hist(E2,xlab="residuals", main="")
#  #OKish
#plot(reach$arabinose.nrr, E2, xlab="NRR", ylab="Residuals")
#  #about the same


#use stream as a random factor
M3<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream, 
        data=reach, na.action=na.omit, method="REML")
anova(M3)
anova(M1,M3)
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

M4<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit)

M5<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit)

M6<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit)
#won't converge even with control=ctrl

M7<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit)

M8<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit)
#won't converge even with control=ctrl

M9<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit)

M10<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf7, data=reach, na.action=na.omit)

anova(M3,M4,M5,M7,M9,M10)
anova(M3,M5)
anova(M5)
summary(M5)
  #M5 is best, but no longer significant

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M11<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf8, data=reach, na.action=na.omit)

M12<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit)

M13<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf10, data=reach, na.action=na.omit)
  #false convergence

anova(M3,M5,M11,M12)
  #M12 is best
anova(M12)

E12<-(residuals(M12))
qqnorm(E12)
ad.test(E12)
  #not normal
hist(E12, xlab="residuals", main="")
  #outliers
plot(reach$arabinose.nrr, E12, xlab="NRR", ylab="Residuals")
  #linear

#try a different alternate variance structure with fitted values
vf11 = varPower(form = ~ fitted(.)|arabinose.nrr)
vf12 = varPower(form = ~ fitted(.)|CBOM.DM)
vf13 = varPower(form = ~ fitted(.)|FBOM.DM)

M14<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf11, data=reach, na.action=na.omit)
M15<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit)
M16<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit)

anova(M5,M12,M14,M15,M16)
  #M12 is still best
  #M14,M15,M16 all have identical LL output...seems fishy

#try combining best variance structures from above (vf2=fitted, and vf9=stream )
M17<-lme(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=varComb(vf2,vf9), data=reach, na.action=na.omit)
anova(M5,M12,M17)
  #M12 still the best with stream as random factor and variance weights for stream

E12<-(residuals(M12))
qqnorm(E12)
ad.test(E12)
  #not normal
hist(E12, xlab="residuals", main="")
  #yuck
plot(reach$arabinose.nrr, E12, xlab="NRR", ylab="Residuals")
  #residuals are still linear or a power function
anova(M12)


#let's try normalizing the data in case that will pull the outlier down into the cloud
reach$L.arabinose.nrr<-log(reach$arabinose.nrr+1)

#start basic
M18<-gls(L.arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM,
         data=reach, na.action=na.omit, method="REML")
#add stream as random
M19<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         data=reach, na.action=na.omit)
#try weight from the best model above
M20<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit)
anova(M18,M19,M20)
  #M20 is best    

E20<-(residuals(M20))
qqnorm(E20)
qqline(E20)
ad.test(E20)  
  #residuals are normal
hist(E20, xlab="residuals", main="")
  #looks OK
plot(reach$L.arabinose.nrr, E20, xlab="NRR", ylab="Residuals")
  #is this better than the unnormalized?  still has 4 large points, but the vertical spread looks better

#refit best model with ML
M21<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit,
         method="ML")
M22<-lme(L.arabinose.nrr~CBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit,
         method="ML")
M23<-lme(L.arabinose.nrr~CBOM.DM+FBOM.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit,
         method="ML")
M24<-lme(L.arabinose.nrr~FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit,
         method="ML")

anova(M21,M22,M23,M24)
#full model is best

#alternatively, drop significant terms from the full model until all terms are significant
summary(M20)
#drop FBOM as least significant
M25<-lme(L.arabinose.nrr~CBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit,
         method="REML")

summary(M25)
anova(M25)
#all terms are significant, so stop here.  same interpretation as using ML:the full model M20
  #containing CBOM and PERI as signifcant but FBOM as not compared to the term deletion M25
  #which deletes FBOM but retains CBOM and PERI as significant

E25<-(residuals(M25))
qqnorm(E25)
qqline(E25)
ad.test(E25)  
  #residuals are normal and better
hist(E25, xlab="residuals", main="")
  #looks OK, but not great
plot(reach$L.arabinose.nrr, E25, xlab="NRR", ylab="Residuals")
  #is this better than the unnormalized?  the vertical spread looks just as bad.
    #compare M12, M20, M25

#interpretation is that the biofilm response to arabinose is higher with greater standing stocks
  #of CBOM but lower with standing stock of periphyton, but the effect size looks very small

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
  #same power function pattern we've been seeing

#optimize random effects before fixed effects
##drop peri.dm from the model since it's not significant and compare the two
#M2<-gls(cellobiose.nrr~CBOM.DM+FBOM.DM,data=reach,na.action=na.omit)
#anova(M2)
#
#anova(M1,M2)
##M2 is better, and FBOM can be dropped
#M3<-gls(cellobiose.nrr~CBOM.DM,data=reach,na.action=na.omit)
#anova(M3)
#
#anova(M2,M3)
##M3 is better

#E3<-(residuals(M3))
#qqnorm(E3)
#qqline(E3)
#ad.test(E3)  
#  #residuals are normal
#hist(E3, xlab="residuals", main="")
#  #OKish
#plot(reach$cellobiose.nrr, E3, xlab="NRR", ylab="Residuals")
#  #exponential relationship

#use stream as a random factor
M4<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream, 
        data=reach, na.action=na.omit, method="REML")
anova(M1,M4)
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

M5<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf1, data=reach, na.action=na.omit)

M6<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf2, data=reach, na.action=na.omit)
  #won't converge even with control=ctrl

M7<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf3, data=reach, na.action=na.omit)
  #won't converge even with control=ctrl

M8<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf4, data=reach, na.action=na.omit)
  #won't converge even with control=ctrl

M9<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf5, data=reach, na.action=na.omit)
  #won't converge even with control=ctrl

M10<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
        weights=vf6, data=reach, na.action=na.omit)
  #won't converge even with control=ctrl

M11<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf7, data=reach, na.action=na.omit)

anova(M4,M5,M11)
  #M4 is still the best

#try different varIdent structures
vf1 = varIdent(form = ~ 1|reach)#pasted from above for organization sake
vf8 = varIdent(form = ~ 1|season)
vf9 = varIdent(form = ~ 1|stream)
vf10 = varIdent(form= ~ 1|reach*season)

M12<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf8, data=reach, na.action=na.omit)

M13<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach, na.action=na.omit)

M14<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
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
vf14 = varPower(form=~fitted(.)|PERI.DM)

M15<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf11, data=reach, na.action=na.omit)
M16<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf12, data=reach, na.action=na.omit)
M17<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf13, data=reach, na.action=na.omit)
M18<-lme(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf14, data=reach, na.action=na.omit)


anova(M13,M15,M16,M17,M18)
  #M13 is best, same strange identical results for vf11-vf14

#let's try normalizing the data
reach$L.cellobiose.nrr<-log(reach$cellobiose.nrr+1)

#start simple
M19<-gls(L.cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, 
         data=reach, na.action=na.omit)
#add stream as random 
M20<-lme(L.cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         data=reach, na.action=na.omit)
#try the weights from the best model above
M21<-lme(L.cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit)
anova(M19,M20,M21)
#M21 has lower AIC     

E21<-(residuals(M21))
qqnorm(E21)
qqline(E21)
ad.test(E21)  
  #residuals are normal
hist(E21, xlab="residuals", main="")
  #looks OK
plot(reach$L.cellobiose.nrr, E21, xlab="NRR", ylab="Residuals")
  #residuals are still messed up, but probably have better vertical spread

anova(M21)
  #significant results
summary(M21)
  #drop CBOM as least significant

M22<-lme(L.cellobiose.nrr~FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit)
summary(M22)
  #drop FBOM
M23<-lme(L.cellobiose.nrr~PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit)
summary(M23)
  #periphyton is significant, but lower stocks are associated with higher biofilm response to cellobiose

#approach using ML
M24<-lme(L.cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit, method="ML")
M25<-lme(L.cellobiose.nrr~CBOM.DM+PERI.DM, random=~1|stream,
         weights=vf9, data=reach,na.action=na.omit, method="ML")
anova(M24,M25)
  #no significant differences with loss of FBOM, but marginaly lower AIC with full model?


#interpretation is that the biofilm response to cellobiose is stronger with lower 
  #standing stocks of CBOM, but the model residuals are screwed up so maybe not valid?
#try to average all carbon sources together because the NDS data showed no carbon effect
  #or do a MANOVA approach to analyze the carbon source simultaneously

##############################################################################
##############################################################################
#multivariate approach
#this might help with model selection  http://www.inside-r.org/packages/cran/MAINT.Data/docs/MANOVA

y <- cbind(reach$glucose.nrr,reach$arabinose.nrr,reach$cellobiose.nrr)
M1 <- manova(y ~ CBOM.DM+FBOM.DM+PERI.DM, data=reach, na.action=na.omit)
summary.aov(M1, test="Pillai")#used to get the univariate statistics
summary(M1,test="Pillai")#used for overall model statistics

residuals(M1)#can I extract each column into a separate vector to look at fits?
M1.r<-data.frame(residuals(M1))
names(M1.r)

fitted(M1)
M1.f<-data.frame(fitted(M1))
names(M1.f)

plot(M1.r$X1,M1.f$X1, xlab="Residuals", ylab="Fitted")
  #not awful
plot(M1.r$X2,M1.f$X2, xlab="Residuals", ylab="Fitted")
  #not awful
plot(M1.r$X3,M1.f$X3, xlab="Residuals", ylab="Fitted")
  #these aren't great

extractAIC(M1)#first value is effective df, second is AIC

#PERI.DM is not significant in the univariate models or the overal model, so drop and rerun
M2 <- manova(y ~ CBOM.DM+FBOM.DM, data=reach, na.action=na.omit)
summary.aov(M2, test="Pillai")#used to get the univariate statistics
summary(M2,test="Pillai")#used for overall model statistics

residuals(M2)#can I extract each column into a separate vector to look at fits?
M2.r<-data.frame(residuals(M2))
names(M2.r)

fitted(M2)
M2.f<-data.frame(fitted(M2))
names(M2.f)

plot(M2.r$X1,M2.f$X1, xlab="Residuals", ylab="Fitted")
#not awful
plot(M2.r$X2,M2.f$X2, xlab="Residuals", ylab="Fitted")
#not awful
plot(M2.r$X3,M2.f$X3, xlab="Residuals", ylab="Fitted")
#these aren't great

extractAIC(M1)
extractAIC(M2)#first value is effective df, second is AIC
#second model has lower AIC

#FBOM.DM is not significant in the overall model, so remove it and run again

M3 <- manova(y ~ CBOM.DM, data=reach, na.action=na.omit)
summary.aov(M3, test="Pillai")#used to get the univariate statistics
summary(M3,test="Pillai")#used for overall model statistics

residuals(M3)#can I extract each column into a separate vector to look at fits?
M3.r<-data.frame(residuals(M3))
names(M3.r)

fitted(M3)
M3.f<-data.frame(fitted(M3))
names(M3.f)

plot(M3.r$X1,M3.f$X1, xlab="Residuals", ylab="Fitted")
#not awful
plot(M3.r$X2,M3.f$X2, xlab="Residuals", ylab="Fitted")
#not awful
plot(M3.r$X3,M3.f$X3, xlab="Residuals", ylab="Fitted")
#these aren't great

extractAIC(M1)#first value is effective df, second is AIC
extractAIC(M2)
extractAIC(M3)
  #second model has lowest AIC
  #final model has CBOM as significant in all univariate summaries and the overall summary

####################################################################
####################################################################

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

#N acquistion enzymes NACE.DM

M1<-gls(NACE.DM~season+reach+stream, data=reach,na.action=na.omit)
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
        weights=vf1, data=reach,na.action=na.omit)
M4<-gls(NACE.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)
M5<-gls(NACE.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(NACE.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)
M7<-gls(NACE.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M6,M7)#accounting for stream variance is best, which makes sense
  #M3 is best


qqnorm(residuals(M3))
qqline(residuals(M3))
ad.test(residuals(M3))  
  #residuals are normal
plot(M3)
  #not bad, but there's an odd gap between 400-900
hist(residuals(M3))
plot(reach$NACE.DM,residuals(M3))
  #looks OK for the n

plot(reach$season,residuals(M3), xlab="Season", ylab="Residuals")
#OK
plot(reach$reach,residuals(M3), xlab="Reach", ylab="Residuals")
#OK
plot(reach$stream,residuals(M3), xlab="Stream", ylab="Residuals")
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
        weights=varComb(vf1,vf7), data=reach,na.action=na.omit)

M10<-gls(NACE.DM~season+reach+stream, 
        weights=varComb(vf1,vf8), data=reach,na.action=na.omit)

M11<-gls(NACE.DM~season+reach+stream, 
        weights=varComb(vf1,vf9), data=reach,na.action=na.omit)

anova(M3,M9,M10,M11)
  #model 3 is still the best 

anova(M3)

#term deletion to optimize model, remove stream which is least significant
  #rerun M3 with ML to compare different fixed effects
M3<-gls(NACE.DM~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit, metho="ML")

M12<-gls(NACE.DM~season+reach, 
        weights=vf1, data=reach, na.action=na.omit, method="ML")

anova(M3,M12)
  #model 12 is better than M3 because no significant difference
anova(M12)

M13<-gls(NACE.DM~season, weights=vf1, 
         data=reach, na.action=na.omit, method="ML")

anova(M12,M13)
  #model 13 is the best because of no significant difference and same AIC

M14<-gls(NACE.DM~season, weights=vf1, 
         data=reach, na.action=na.omit, method="REML")

anova(M14)
summary(M14)

qqnorm(residuals(M14))
qqline(residuals(M14))
ad.test(residuals(M14))  
  #residuals are not normal
plot(M14)
hist(residuals(M14))
plot(reach$NACE.DM,residuals(M14))
  #looks OK for the n

plot(reach$season,residuals(M14), xlab="Season", ylab="Residuals")
  #OK
plot(reach$reach,residuals(M14), xlab="Reach", ylab="Residuals")
  #OK
plot(reach$stream,residuals(M14), xlab="Stream", ylab="Residuals")
  #OK

  #interpretation is that Fall has the most N acquiring enzymes NACE compared to spring and summer
    #possibly due to a pulse of recalcitrant terrestrial C?  
  #M14 has somewhat worse looking diagnostics, but we deleted the insignificant terms 
    #and the interpretation is the same

####################################################################
####################################################################

#N acquistion enzymes NACE.C

M1<-gls(NACE.C~season+reach+stream, data=reach,na.action=na.omit)
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
plot(reach$NACE.C,residuals(M1))
  #positive relationship?
y<-lm(residuals(M1)~reach$NACE.C)
summary(y)
  #yep

plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #OK
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #OK, but daylight has bigger spread
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK, stream looks better than with NACE.DM

#looks like stream is weirdest. let's try it as a random factor
M1<-gls(NACE.C~season+reach+stream, 
        data=reach,na.action=na.omit, method="ML")#rerun with ML for comparison of different fixed effects
M2<-lme(NACE.C~season+reach, random=~1|stream, 
        data=reach,na.action=na.omit, method="ML")


anova(M1,M2)
#M1 is better, so maybe nesting is not the best here

M1<-gls(NACE.C~season+reach+stream, data=reach,na.action=na.omit, method="REML")

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varIdent(form=~1|stream*season*reach)#this one won't work because there are 0 df for the three interaction

M3<-gls(NACE.C~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit)
M4<-gls(NACE.C~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)
M5<-gls(NACE.C~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(NACE.C~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)
M7<-gls(NACE.C~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
  #false convergence

anova(M1,M3,M4,M5,M6)#accounting for stream variance is best, which makes sense
  #M3 is best


qqnorm(residuals(M3))
qqline(residuals(M3))
ad.test(residuals(M3))  
  #residuals are normal
plot(M3)
  #not bad
hist(residuals(M3))
  #looks better
plot(reach$NACE.C,residuals(M3))
  #looks OK for the n

plot(reach$season,residuals(M3), xlab="Season", ylab="Residuals")
  #little weirder
plot(reach$reach,residuals(M3), xlab="Reach", ylab="Residuals")
  #OK
plot(reach$stream,residuals(M3), xlab="Stream", ylab="Residuals")
  #not as good as M1

#by stream is still funky so try stream as a random effect along with the varIdent
M8<-lme(NACE.C~season+reach+stream, random=~1|stream,
        weights=vf1, data=reach,na.action=na.omit)

anova(M3,M8)
  #M3 still wins

#try a model based on fitted values
vf7 = varPower(form = ~ fitted(.))
vf8 = varExp(form = ~ fitted(.))
vf9 = varConstPower(form = ~ fitted(.))

M9<-gls(NACE.C~season+reach+stream, 
        weights=varComb(vf1,vf7), data=reach,na.action=na.omit)

M10<-gls(NACE.C~season+reach+stream, 
         weights=varComb(vf1,vf8), data=reach,na.action=na.omit)

M11<-gls(NACE.C~season+reach+stream, 
         weights=varComb(vf1,vf9), data=reach,na.action=na.omit)

anova(M3,M9,M10,M11)
  #model 3 is still the best 

anova(M3)

#term deletion to optimize model, remove stream which is least significant
#rerun M3 with ML to compare different fixed effects
M3<-gls(NACE.C~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit, method="ML")

M12<-gls(NACE.C~season+stream, 
         weights=vf1, data=reach, na.action=na.omit, method="ML")

anova(M3,M12)
#model 12 is better than M3 because no significant difference
anova(M12)

M13<-gls(NACE.C~season, weights=vf1, 
         data=reach, na.action=na.omit, method="ML")

anova(M12,M13)
#model 13 is the best because of no significant difference and better AIC

M14<-gls(NACE.C~season, weights=vf1, 
         data=reach, na.action=na.omit, method="REML")

anova(M14)
summary(M14)

qqnorm(residuals(M14))
qqline(residuals(M14))
ad.test(residuals(M14))  
  #residuals are not normal
plot(M14)
  #horizontal spread is categorical, but vertical spread is OK
hist(residuals(M14))
  #not as good
plot(reach$NACE.C,residuals(M14))
  #crazy!

plot(reach$season,residuals(M14), xlab="Season", ylab="Residuals")
  #OK
plot(reach$reach,residuals(M14), xlab="Reach", ylab="Residuals")
  #OK
plot(reach$stream,residuals(M14), xlab="Stream", ylab="Residuals")
  #this isn't an improvement

#interpretation is that Fall has the most N acquiring enzymes NACE compared to spring and summer
#possibly due to a pulse of recalcitrant terrestrial C?  
#M14 has somewhat worse looking diagnostics, but we deleted the insignificant terms 
#and the interpretation is the same

#same interpretation as with NACE.DM, but with same diagnostic problems.  

#stick with the model with deleted terms or more complicated model with better diagnostics?  
  #the interpretation is the same either way.

####################################################################
####################################################################

#N acquisition -- LEU.C produces models with no significant factors

M1<-gls(LEU.DM~season+reach+stream, data=reach,na.action=na.omit, method="REML")
anova(M1)
summary(M1)
  #almost higher in spring compared to fall p=0.0503

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #not normal
plot(M1)
  #variation proportional to fit
hist(residuals(M1))
plot(reach$LEU.DM,residuals(M1))
  #positive relationship?
  y<-lm(residuals(M1)~reach$LEU.DM)
  summary(y)
    #by outlier...need to deal with this

#try to use stream as a random factor
M1<-gls(LEU.DM~season+reach+stream, data=reach,na.action=na.omit, method="ML")
M2<-lme(LEU.DM~season+reach, random=~1|stream, 
          data=reach,na.action=na.omit, method="ML")
anova(M1,M2)
  #M1 is better, so the random factor is not the best

M1<-gls(LEU.DM~season+reach+stream, data=reach,na.action=na.omit, method="REML")

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

ctrl<-glsControl(returnObject=TRUE, maxIter=5200, msMaxIter=5200)

M3<-gls(LEU.DM~season+reach+stream, 
       weights=vf1, data=reach,na.action=na.omit)  
M4<-gls(LEU.DM~season+reach+stream, 
        weights=vf2, data=reach,na.action=na.omit)  
M5<-gls(LEU.DM~season+reach+stream, 
        weights=vf3, data=reach,na.action=na.omit)
M6<-gls(LEU.DM~season+reach+stream, 
        weights=vf4, data=reach,na.action=na.omit)
M7<-gls(LEU.DM~season+reach+stream, 
        weights=vf5, data=reach,na.action=na.omit)
M8<-gls(LEU.DM~season+reach+stream, 
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(LEU.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)
  #no convergence

anova(M1,M3,M4,M5,M6,M7,M8)
  #M8 is best, using power of fitted as an alternate variance structure

#try a combined variance structure just for fun
M10<-gls(LEU.DM~season+reach+stream, 
        weights=varComb(vf6,vf4), data=reach,na.action=na.omit)  
  #no convergence

anova(M8,M10)
  #M8 wins based on AIC

anova(M5)
summary(M5)
  #no significant factors except:
    #M5 with reach and stream
    #M6 and M7 has freaky F values
  #shows daylight having more N acquiring enzymes than buried
    #given how flaky this analysis is, I'm not inclined to include it

plot(reach$season,residuals(M5), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M5), reach$season)  
  #not OK
plot(reach$reach,residuals(M5), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M5), reach$reach)  
  #not OK
plot(reach$stream,residuals(M5), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M5), reach$stream)  
  #not OK

qqnorm(residuals(M5))
qqline(residuals(M5))
ad.test(residuals(M5))  
  #not OK
plot(M5)
  #these look better
hist(residuals(M5))
plot(reach$LEU.DM,residuals(M5))
  #positive relationship?
  y<-lm(residuals(M5)~reach$LEU.DM)
  summary(y)
    #driven by the same outlier

#This analysis looks like a marginal analysis at best because the best model on AIC
  #has no significant factors and the model with interesting significant factors
  #doesn't pass diagnostics

##################################################################
##################################################################

#Sulfur acquisition  (not going to focus on this analysis in favor of the C,N story)

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
  #buried has big spread
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
  #OK, but Este looks tight compared to Amberly

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
        weights=vf6, data=reach,na.action=na.omit)
  #no convergence
M9<-gls(PHOS.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)
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
  #unreasonable F values.  need to look at other models

anova(M10)
  #M3 has season and stream as significant
  #M4, M5 have no interesting significant factors
  #M6, M7, M10 have unreasonable F values

#let's see if graphical validation looks OK for M3

plot(reach$season,residuals(M3), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M3), reach$season)  
  #OK
plot(reach$reach,residuals(M3), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M3), reach$reach)  
  #OK
plot(reach$stream,residuals(M3), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M3), reach$stream)  
  #not great
  
qqnorm(residuals(M3))
qqline(residuals(M3))
ad.test(residuals(M3))  
  #not good, but it could be worse
plot(M3)
  #this looks OK
hist(residuals(M3))
plot(reach$PHOS.DM,residuals(M3))
  #looks better when we started, but those two outliers are still odd

anova(M3)
summary(M3)
#if we use any model from this, it should be either M1 or M3.  M3 has better
#graphical diagnostics, and just marginally worse AIC. Interpretation is that
#spring and summer have lower P acquistion enzymes compared to fall

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
        weights=vf6, data=reach,na.action=na.omit)
  #false convergence
M9<-gls(DOPAH2.DM~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M6,M9)
  #M1 is best on AIC, none are significantly different
  #graphical analysis indicates M4 or M5 might handle residual variation better 

#try a combined variance structure just for fun
M10<-gls(DOPAH2.DM~season+reach+stream, 
         weights=varComb(vf2,vf7), data=reach,na.action=na.omit)  
anova(M1,M4,M5,M10)

anova(M10)
  #M3, M4, M5, M9 reach is significant
  #M6 has unreasonable F values
  #M10 has nothing significant

plot(reach$season,residuals(M9), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M9), reach$season)  
  #OK
plot(reach$reach,residuals(M9), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M9), reach$reach)  
  #OK
plot(reach$stream,residuals(M9), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M9), reach$stream)  
  #not great
  
qqnorm(residuals(M9))
qqline(residuals(M9))
ad.test(residuals(M9))  
  #not good
plot(M9)
  #this looks OK
hist(residuals(M9))
plot(reach$DOPAH2.DM,residuals(M9))

  #M5 does the overall best job
anova(M5)
summary(M5)
  #interpretation is that daylight has lower effort to acquire recalcitrant carbon compared to 
    #buried

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
        weights=vf6, data=reach,na.action=na.omit)
  #no convergence
M9<-gls(POX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M6,M9)
  #M5 has lowest AIC with variance accounting for season

#try a combined variance structure just for fun
M10<-gls(POX~season+reach+stream, 
         weights=varComb(vf3,vf7), data=reach, na.action=na.omit) 

anova(M1,M5,M10)
  #the combined model is not an improvement

qqnorm(residuals(M5))
qqline(residuals(M5))
ad.test(residuals(M5))  
  #not so good

plot(M5)
  #not horrible
hist(residuals(M5))
  #major outliers
plot(reach$POX,residuals(M5))
  #outliers

plot(reach$season,residuals(M5), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M5), reach$season)  
  #not OK big variation in spring
plot(reach$reach,residuals(M5), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M5), reach$reach)  
  #OK
plot(reach$stream,residuals(M5), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M5), reach$stream)  
  #OK

#try normalizing?
reach$L.POX<-log(reach$POX+1)

#start simple
M11<-gls(L.POX~season+reach+stream, data=reach,na.action=na.omit)
anova(M11)
summary(M11)

#try best alternate variance structure from above
M12<-gls(L.POX~season+reach+stream, data=reach, na.action=na.omit,
         weights=vf3)

anova(M11,M12)
#M11 looks like the best best

qqnorm(residuals(M11))
qqline(residuals(M11))
ad.test(residuals(M11))  
  #good

plot(M11)
  #much better
hist(residuals(M11))
  #still an outlier, but better
plot(reach$POX,residuals(M11))
  #outliers still there but better horizontal alignment

plot(reach$season,residuals(M11), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M11), reach$season)  
  #OK spring variation looks better
plot(reach$reach,residuals(M11), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M11), reach$reach)  
  #OK
plot(reach$stream,residuals(M11), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M11), reach$stream)  
  #OK

#graphical analysis shows a big improvement in the normalized data

anova(M11)

#term deletion 
M11.full<-gls(L.POX~season+reach+stream, data=reach,na.action=na.omit,
         method="ML")
anova(M11.full)

M12<-gls(L.POX~season+reach, data=reach,na.action=na.omit,
           method="ML")
anova(M11.full,M12)
  #deleting stream is OK based on n.s. p value and lower AIC

anova(M12)
M13<-gls(L.POX~reach, data=reach,na.action=na.omit,
         method="ML")
anova(M12,M13)
  #deleting season is OK based on n.s. p value and lower AIC

M13.full<-gls(L.POX~reach, data=reach,na.action=na.omit,
              method="REML")
anova(M13.full)
summary(M13.full)

  #daylighted reaches have lower effort invested in acquiring recalcitrant carbon

##################################################################
##################################################################
#LCI is an index of carbon lability calculcated as recalcitrant/(labile+recalcitrant)
  #so a smaller number means there is more labile carbon and a bigger number means more
  #recalcitrant carbon

M1<-gls(LCI~season+reach+stream, data=reach, na.action=na.omit)
anova(M1)
summary(M1)
  #smaller number in daylight compared to buried (more labile in daylight),
    #but the seasonal pattern is opposite what we might expect.
    #this is a ratio rather than a raw number

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #OK

plot(M1)
  #unimodal or not bad?
hist(residuals(M1))
plot(reach$LCI,residuals(M1))
  #increasing variation with LCI

plot(reach$season,residuals(M1), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M1), reach$season)  
  #looks great
plot(reach$reach,residuals(M1), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M1), reach$reach)  
  #looks great
plot(reach$stream,residuals(M1), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M1), reach$stream)  
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
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(LCI~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M6,M8,M9)
  #M6 which is stream*reach is best
anova(M6)
  #F values are unrealistic

anova(M8)
  #M3, M8, M9 reach is significant
  #M4, M5, season, reach, stream significant

qqnorm(residuals(M5))
qqline(residuals(M5))
ad.test(residuals(M5))  
#OK

plot(M10)
  #better
hist(residuals(M10))
plot(reach$LCI,residuals(M10))
  #increasing variation with LCI

plot(reach$season,residuals(M10), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M10), reach$season)  
#looks great
plot(reach$reach,residuals(M10), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M10), reach$reach)  
#looks great
plot(reach$stream,residuals(M10), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M10), reach$stream)  
#OK

anova(M4)
summary(M4)
  #M4 seems to handle the residuals versus observations the best.  
  #daylight has lower values (greater carbon lability) than buried, but spring/summer have
    #higher values (lower carbon lability) than  fall?

#try a combined variance structure just for fun
M10<-gls(LCI~season+reach+stream, 
         weights=varComb(vf2,vf6), data=reach, na.action=na.omit) 

anova(M4,M10)
anova(M10)

#the combined structure doesn't handle the trend of increasing residuals with observations as well, 
  #but the model output is more intuitive for interpreting the seasonal effect

##################################################################
##################################################################
#HIX is from Pennino's work, it is the humification index, so a bigger number
  #means more humic compounds

M1<-gls(HIX~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #almost lower in spring compared to fall, p=0.0538

qqnorm(residuals(M13))
qqline(residuals(M13))
  #my goodness
ad.test(residuals(M13))  
  #how is this not worse than p=0.03?

plot(M13)
  #one major outlier that is underestimated
hist(residuals(M13))
plot(new.reach$HIX,residuals(M13))
  #one distinct outlier

plot(new.reach$season,residuals(M13), xlab="Season", ylab="Residuals")
  #there's that outlier
bartlett.test(residuals(M13), new.reach$season)  
  #OK
plot(new.reach$reach,residuals(M13), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M13), new.reach$reach)  
  #not OK
plot(new.reach$stream,residuals(M13), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M13), new.reach$stream)  
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
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(HIX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M6,M7,M8,M9)
  #M3 wins, stream
anova(M3)
  #graphical diagnostics don't show an improvement with the outlier.  we might have to dump it?

#try normalizing
reach$L.HIX<-log(reach$HIX)

#start simple
M10<-gls(L.HIX~season+reach+stream, data=reach,na.action=na.omit)

#try best alternate variance structure from above
M11<-gls(L.HIX~season+reach+stream, data=reach, na.action=na.omit,
         weights=vf1)

anova(M10,M11)
  #M11 wins

#sub M11 into the plot functions above...no fix on the outlier.  Try removing it
new.reach<-reach[reach$HIX>0.71,]

M12<-gls(HIX~season+reach+stream, data=new.reach, na.action=na.omit)
anova(M1)
summary(M1)

M13<-gls(HIX~season+reach+stream, 
        weights=vf1, data=new.reach,na.action=na.omit)  
M14<-gls(HIX~season+reach+stream, 
        weights=vf2, data=new.reach,na.action=na.omit)  
M15<-gls(HIX~season+reach+stream, 
        weights=vf3, data=new.reach,na.action=na.omit)
M16<-gls(HIX~season+reach+stream, 
        weights=vf4, data=new.reach,na.action=na.omit)  
M17<-gls(HIX~season+reach+stream, 
        weights=vf5, data=new.reach,na.action=na.omit)
  #false convergence
M18<-gls(HIX~season+reach+stream, 
        weights=vf6, data=new.reach,na.action=na.omit)
M19<-gls(HIX~season+reach+stream, 
        weights=vf7, data=new.reach,na.action=na.omit)

anova(M12,M13,M14,M15,M16,M18,M19)
  #M13 wins

anova(M13)
  #sub M13 into graphical diagnostics above
    #outlier removed in new.reach looks a lot better, but stream is still a little goofy

summary(M13)
  #interpretation:  spring and summer have lower humification than fall, possibly due to less 
    #terrestrial input.  daylight has higher humification than buried, possibly because 
    #terrestrial inputs occur in daylighted reaches but not buried reaches?
  


##################################################################
##################################################################
#BIX is from Pennino's work, it is the biological freshness index.  Bigger number
#means more autochthony

M1<-gls(BIX~season+reach+stream, data=reach,na.action=na.omit)
anova(M1)
summary(M1)
  #spring and summer have higher BIX than fall

qqnorm(residuals(M4))
qqline(residuals(M4))
ad.test(residuals(M4))  
  #good
hist(residuals(M4))
plot(reach$BIX,residuals(M4))
  #major categorical variation, but vertical spread not awful

plot(reach$season,residuals(M4), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M4), reach$season)  
  #OK
plot(reach$reach,residuals(M4), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M4), reach$reach)  
  #OK
plot(reach$stream,residuals(M4), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M4), reach$stream)  
  #OK

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

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
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(BIX~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M8,M9)
anova(M9)
  #M3 has season and reach as significant
  #M4 has season as significant  
  #M5, M8, M9 have unrealistic F values

  #M8 wins, but graphical analysis shows it doesn't handle season well
  #M3 probably handles the residuals versus observations a bit better
  #is this as good as it gets given the bimodality?  I can't imagine
    #log normalization will help given the bimodality
  
anova(M3)
summary(M3)

#term deletion
M3.full<-gls(BIX~season+reach+stream, 
        weights=vf1, data=reach,na.action=na.omit,
        method="ML")  

M10<-gls(BIX~season+reach, 
         weights=vf1, data=reach,na.action=na.omit,
         method="ML")  

anova(M3.full, M10)
  #M10 has lower AIC but not different, so it's better

anova(M10)
  #all terms are significant, so we're good.  Rerun with REML 
M10.full<-gls(BIX~season+reach, 
         weights=vf1, data=reach,na.action=na.omit,
         method="REML")  
anova(M10.full)
summary(M10.full)

  #spring and summer have a higher freshness index than fall
  #daylight has less freshness than than buried (possibly due to green/litter fall inputs?)
  #but effect size is small and barely significant

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
        weights=vf6, data=reach,na.action=na.omit)
M9<-gls(FI~season+reach+stream, 
        weights=vf7, data=reach,na.action=na.omit)

anova(M1,M3,M4,M5,M7,M8,M9)
  #looks like M9 wins.

anova(M9)
  #M3, M4, M5, M8, M9 season is significant
  #M7 has unreasonable F values
  #sub these models into the graphical evaluation above
    #looks like M3-5 might handle the data the best, but the bimodal nature
      #is a challenge
anova(M5)

summary(M5)
  #spring and summer have a higher FI than fall indicating more microbial compared to terrestrial
  

  #random thought...is high fall variance compared to narrow summer variance
    #accounted for by different disturbance frq (more t-storm summer) or
    #differences in OM (leaf input in fall)? Maybe we need a covariate like CBOM standing stock?

  #need to decide which model is best and then do backward selection

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

M1<-gls(P2H~season+reach+stream, data=new.reach,na.action=na.omit)
anova(M1)
summary(M1)  
  #higher in spring and summer compared to fall

qqnorm(residuals(M9))
qqline(residuals(M9))
ad.test(residuals(M9))  
  #that one outlier is bad, but otherwise OK
hist(residuals(M9))
plot(new.reach$P2H,residuals(M9))
  #bimodal spread with more variation in summer/spring

plot(new.reach$season,residuals(M9), xlab="Season", ylab="Residuals")
bartlett.test(residuals(M9), new.reach$season)  
  #not good
plot(new.reach$reach,residuals(M9), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(M9), new.reach$reach)  
  #OK
plot(new.reach$stream,residuals(M9), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(M9), new.reach$stream)  
  #spring buried eastgate, which was the same bit outlier earlier

#try alternate variance structures
vf1<-varIdent(form = ~1|stream)
vf2<-varIdent(form=~1|reach)
vf3<-varIdent(form=~1|season)
vf4<-varIdent(form=~1|stream*reach)
vf5<-varIdent(form=~1|stream*season)
vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

M3<-gls(P2H~season+reach+stream, 
        weights=vf1, data=new.reach,na.action=na.omit)  
M4<-gls(P2H~season+reach+stream, 
        weights=vf2, data=new.reach,na.action=na.omit)  
M5<-gls(P2H~season+reach+stream, 
        weights=vf3, data=new.reach,na.action=na.omit)
M6<-gls(P2H~season+reach+stream, 
        weights=vf4, data=new.reach,na.action=na.omit)  
  #false convergence
M7<-gls(P2H~season+reach+stream, 
        weights=vf5, data=new.reach,na.action=na.omit)
M8<-gls(P2H~season+reach+stream, 
        weights=vf6, data=new.reach,na.action=na.omit)
M9<-gls(P2H~season+reach+stream, 
        weights=vf7, data=new.reach,na.action=na.omit)

anova(M1,M3,M4,M5,M7,M8,M9)
  #M7 look like it wins, stream*season

anova(M9)
  #M4 has season significant
  #M3, M5, M9 has season and stream significant
  #M8 has everything significant
  #M7 has unreasonable F ratios

#that's the same outlier we did the reduced data set for rerun with data=new.reach

anova(M1,M3,M4,M5,M8,M9)
  #looks like M5, M8, or M9 might be the best
  #graphical diagnostics all look about the same

summary(M8)
  #spring and summer have the highest P2H in all models, one model shows daylight with less
    #P2H but a very small effect size

##################################################################
##################################################################
#CQI and LCI should be inversely related to each other
  #bigger CQI=better carbon, bigger LCI=worse carbon

M1<-gls(CQI~LCI, data=new.reach, na.action=na.omit)
anova(M1)
summary(M1)  
  #negative relationship

plot(log(new.reach$CQI),new.reach$LCI)
#gack...need to get rid of that outlier...Eastgate Fall Daylight

#remove outlier
new2.reach<-reach[reach$CQI<20,]

M1<-gls(CQI~LCI, data=new2.reach, na.action=na.omit)
anova(M1)
summary(M1)  
#negative relationship

plot(new2.reach$CQI,new2.reach$LCI)
  #looks better, but non-linear

qqnorm(residuals(M1))
qqline(residuals(M1))
ad.test(residuals(M1))  
  #OK
plot(M1)
  #very distinct pattern
hist(residuals(M1))
  #not good

vf6<-varPower(form = ~ fitted(.))
vf7<-varExp(form = ~ fitted(.))

M2<-gls(CQI~LCI, data=new2.reach, na.action=na.omit,
        weights=vf6)
M3<-gls(CQI~LCI, data=new2.reach, na.action=na.omit,
        weights=vf7)

anova(M1,M2,M3)

#both are significant improvements based on AIC
qqnorm(residuals(M3))
qqline(residuals(M3))
ad.test(residuals(M3))  
#OK
plot(M3)
#very distinct pattern
hist(residuals(M3))
#not good

#those models improved the residuals v. fitted, but not much else
  #do we need to use some non-linear regression technique or just not worry about it?
M4<-nls(new2.reach$CQI~new2.reach$LCI)
anova(M4)

