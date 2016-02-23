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
  #positive relationship with CBOM.DM and FBOM.DM, weaker with arabinose and cellobiose

x<-gls(glucose.nrr~daily.par,data=reach,na.action=na.omit)
  #positive relationship with daily PAR

x<-gls(BIX~season+reach+stream, data=reach,na.action=na.omit)
  #spring and summer have higher BIX than fall

x<-gls(NACE.DM~season+reach+stream, data=reach,na.action=na.omit)
  #significantly lower in spring compared to Fall

x<-gls(LEU.C~season+reach+stream, data=reach,na.action=na.omit)
  #almost higher in spring compared to fall

x<-gls(SULF.DM~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring than fall or summer, same with SULF.C

x<-gls(PHOS.DM~season+reach+stream, data=reach,na.action=na.omit)
  #lower in summer compared to fall

x<-gls(DOPAH2.DM~season+reach+stream, data=reach,na.action=na.omit)
  #lower in daylight compared to buried

x<-gls(DOPA~season+reach+stream, data=reach,na.action=na.omit)
  #almost lower in daylight compared to buried

x<-gls(POX~season+reach+stream, data=reach,na.action=na.omit)
  #almost lower in daylight compared to buried

x<-gls(LCI~season+reach+stream, data=reach, na.action=na.omit)
  #lower in daylight compared to buried

x<-gls(HIX~season+reach+stream, data=reach,na.action=na.omit)
  #almost lower in spring compared to fall
  
x<-gls(BIX~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall

x<-gls(FI~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall

x<-gls(PRO~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall

x<-gls(FUL~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall and in eastgate compared to Amberly

x<-gls(HUM~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall and in eastgate compared to Amberly

x<-gls(P2H~season+reach+stream, data=reach,na.action=na.omit)
  #higher in spring and summer compared to fall

x<-gls(CQI~LCI, data=reach, na.action=na.omit)
  #negative relationship between CQI and LCI

x<-gls(glucose.nrr~P2H, data=reach, na.action=na.omit)
summary(x)

