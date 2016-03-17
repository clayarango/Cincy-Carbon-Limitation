#21-Feb-16
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
library(nlme)
library(nortest)

#attach and evaluate data
reach<-read.table(file="cincy.reach.csv", header=T, sep=',')
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

x<-gls(glucose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit)
summary(x)
x<-gls(arabinose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit)
summary(x)
x<-gls(cellobiose.nrr~CBOM.DM+FBOM.DM+PERI.DM,data=reach,na.action=na.omit)
summary(x)  
  #positive relationship between glucose and CBOM.DM and FBOM.DM, weaker with arabinose, cellobiose p<0.10

x<-gls(glucose.nrr~daily.par,data=reach,na.action=na.omit)
summary(x)  
x<-gls(arabinose.nrr~daily.par,data=reach,na.action=na.omit)
summary(x)  
x<-gls(cellobiose.nrr~daily.par,data=reach,na.action=na.omit)
summary(x)  
  #weak positive relationship between glucose daily PAR, p<0.10 with arabinose, p<0.15 cellobiose


x<-gls(NACE.DM~season+reach+stream, data=reach,na.action=na.omit)#N acq enzymes
summary(x)
  #significantly lower in spring compared to Fall

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
plot(x)
hist(residuals(x))
plot(reach$NACE.DM,residuals(x))
  #positive relationship?
  y<-lm(residuals(x)~reach$NACE.DM)
  summary(y)
  #yep
    
plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK



x<-gls(LEU.C~season+reach+stream, data=reach,na.action=na.omit)
summary(x)  
  #almost higher in spring compared to fall p=0.0503

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #OK, but they don't look great
plot(x)
  #variation proportional to fit
hist(residuals(x))
plot(reach$LEU.DM,residuals(x))
  #positive relationship?
  y<-lm(residuals(x)~reach$NACE.DM)
  summary(y)
    #almost...driven by the same outlier. probably need to deal with this

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #not OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #OK, but they don't look great
plot(x)
  #variation proportional to fit
hist(residuals(x))
plot(reach$LEU.DM,residuals(x))
  #positive relationship?
  y<-lm(residuals(x)~reach$NACE.DM)
  summary(y)
    #almost...driven by the same outlier. probably need to deal with this

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #not OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK



x<-gls(SULF.DM~season+reach+stream, data=reach,na.action=na.omit)
summary(x)
  #higher in spring than fall or summer, same with SULF.C

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #OK, but they don't look great
plot(x)
  #gack
hist(residuals(x))
plot(reach$SULF.DM,residuals(x))
  #three major outliers

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK, but those three big outliers might be from Amberly



x<-gls(PHOS.DM~season+reach+stream, data=reach,na.action=na.omit)
summary(x)
  #lower in summer compared to fall

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #OK, but they don't look great
plot(x)
  #OK
hist(residuals(x))
plot(reach$PHOS.DM,residuals(x))
  #two outliers or increasing trend?

  y<-lm(residuals(x)~reach$PHOS.DM)
  summary(y)
    #significant, but driven by outliers

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK, but ugly
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK, outliers in Amberly again?



x<-gls(DOPAH2.DM~season+reach+stream, data=reach,na.action=na.omit)
summary(x)
  #lower in daylight compared to buried

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #not so good
plot(x)
  #sort of a U shape
hist(residuals(x))
  #major outliers
plot(reach$DOPAH2.DM,residuals(x))
  #increasing trend?

  y<-lm(residuals(x)~reach$DOPAH2.DM)
  summary(y)
    #significant

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK, but ugly
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK?  Este is all over the map, and big outliers in the others



x<-gls(DOPA~season+reach+stream, data=reach,na.action=na.omit)
summary(x)  
  #almost lower in daylight compared to buried, p<0.10
  #we already have DOPAH2, so maybe won't pursue this since it's on the edge anyway

x<-gls(POX~season+reach+stream, data=reach,na.action=na.omit)#phenol oxidase, recalcitrant C
summary(x)
  #almost lower in daylight compared to buried

qqnorm(residuals(x))
  #exponential
qqline(residuals(x))
ad.test(residuals(x))  
  #not so good

plot(x)
  #couple big overestimates
hist(residuals(x))
  #major outliers
plot(reach$POX,residuals(x))
  #increasing trend? or outliers?

  y<-lm(residuals(x)~reach$POX)
  summary(y)
    #significant

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #not OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK

vf<-varIdent(form=~1|season)
x2<-gls(POX~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
  #daylight lower than buried
anova(x,x2)
  #no significant difference with models
  #need more diagnostics here...the model with season variance still has the same problems as without

x<-gls(LCI~season+reach+stream, data=reach, na.action=na.omit)#index of carbon lability (bigger more labile)
summary(x)
  #lower in daylight compared to buried

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #OK

plot(x)
  #unimodal?
hist(residuals(x))
plot(reach$LCI,residuals(x))
  #increasing variation with LCI

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK




x<-gls(HIX~season+reach+stream, data=reach,na.action=na.omit)#humification index
summary(x)
  #almost lower in spring compared to fall, p=0.0538

qqnorm(residuals(x))
qqline(residuals(x))
  #my goodness
ad.test(residuals(x))  
  #how is this not worse than p=0.03?

plot(x)
  #one major outlier that is underestimated
hist(residuals(x))
plot(reach$LCI,residuals(x))
  #increasing variation with LCI and/or downward trend?

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
  #there's that outlier
bartlett.test(residuals(x), reach$season)  
  #OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #not OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK, eastgate, buried spring with the outlier...better check the raw value

vf<-varIdent(form=~1|stream)
x2<-gls(HIX~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2)
  #model with stream is better...need diagnostics



x<-gls(BIX~season+reach+stream, data=reach,na.action=na.omit) #biological freshness of DOM
summary(x)
  #spring and summer have higher BIX than fall
qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
plot(x)
hist(residuals(x))
plot(reach$BIX,residuals(x))
  #major categorical variation

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
#OK
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
#OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
#OK


x<-gls(FI~season+reach+stream, data=reach,na.action=na.omit)#bigger=microbial, lower=terrestrial
summary(x)
  #higher in spring and summer compared to fall

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
plot(x)
hist(residuals(x))
  #big outlier
plot(reach$FI,residuals(x))
  #major categorical variation

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
  #egad
bartlett.test(residuals(x), reach$season)  
  #OK on p value, but ugly
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #OK

vf<-varIdent(form=~1|season)
x2<-gls(FI~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2) #better model accounting for season variance, but need diagnostics 
  #random thought...is high fall variance compared to narrow summer variance
    #accounted for by different disturbance frq (more t-storm summer) or
    #differences in OM (leaf input in fall)? Maybe we need a covariate like days 
    #since last storm or CBOM standing stock?


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


x<-gls(P2H~season+reach+stream, data=reach,na.action=na.omit)
summary(x)  
  #higher in spring and summer compared to fall
  #this should be the better y variable compared to HUM or PRO alone

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #that one outlier is bad, but otherwise OK
plot(x)
  #sort of bimodal, variation increases with fit
hist(residuals(x))
plot(reach$P2H,residuals(x))
  #not awful...just that one outlier

plot(reach$season,residuals(x), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x), reach$season)  
  #not good
plot(reach$reach,residuals(x), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x), reach$reach)  
  #OK
plot(reach$stream,residuals(x), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x), reach$stream)  
  #spring buried eastgate

vf<-varIdent(form=(~1|season*stream))
x2<-gls(P2H~season+reach+stream, data=reach,na.action=na.omit, weights=vf)
summary(x2)
anova(x,x2)
  #much better model with very significant difference between daylight and buried

qqnorm(residuals(x2))
qqline(residuals(x2))
ad.test(residuals(x2))  
  #that outlier is still bad
plot(x2)
  #sort of bimodal, but variation no long increases with fit
hist(residuals(x2))
plot(reach$P2H,residuals(x2))
  #not awful...just that one outlier

plot(reach$season,residuals(x2), xlab="Season", ylab="Residuals")
bartlett.test(residuals(x2), reach$season)  
  #not good
plot(reach$reach,residuals(x2), xlab="Reach", ylab="Residuals")
bartlett.test(residuals(x2), reach$reach)  
  #not good
plot(reach$stream,residuals(x2), xlab="Stream", ylab="Residuals")
bartlett.test(residuals(x2), reach$stream)  
  #that one outlier is still a problem, but everything got tighter



x<-gls(CQI~LCI, data=reach, na.action=na.omit)
summary(x)  
  #negative relationship between CQI and LCI
  #makes sense bc bigger CQI=better carbon, bigger LCI=worse carbon

qqnorm(residuals(x))
qqline(residuals(x))
ad.test(residuals(x))  
  #that one outlier is bad, but otherwise OK
plot(x)
  #very distinct pattern...this looks suspicious
hist(residuals(x))
  #not good
plot(reach$CQI,reach$LCI)
#gack...need to get rid of that outlier...i don't think CQI should be greater than 1

