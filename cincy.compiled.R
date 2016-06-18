#9-May-16
#This R script compiles all our final analyses from cincy.reach.R
#and cincy.nds.R.  Those two files contain the variance optimization 
#and model selection code.

############################################################################################
#load and set up data

#packages
install.packages("nlme")
install.packages("nortest")
install.packages("dplyr")
install.packages("multcomp")
install.packages("MASS")
install.packages("ggplot2")
library(nlme)
library(nortest)
library(dplyr)
library(multcomp)
library(MASS)
library(ggplot2)

#load and evaluate NDS data
nds<-read.table(file="cincy.nds.csv", header=T, sep=',')
str(nds)
unique(nds$stream)
unique(nds$season)
unique(nds$reach)
unique(nds$carbon)

nds$carbon <- factor(nds$carbon, levels = c("Control", "Arabinose", "Cellobiose", "Glucose"))

#load and evaluate REACH data
reach<-read.table(file="cincy.reach.csv", header=T, sep=',')
unique(reach$stream)
unique(reach$season)
unique(reach$reach)
str(reach)

############################################################################################
#Carbon limitation analysis with NDS data
fall<-subset(nds, season=="Fall")
spring<-subset(nds, season=="Spring")
fall.spring<-rbind(fall,spring)
fall.spring<-droplevels(fall.spring)#removes summer as a level bc it has no values

vf1 = varIdent(form = ~1|season*reach) 
vf2 = varExp(form = ~ fitted (.)|season)
vf3 = varIdent(form = ~1|stream)

M.nds<-lme(nrr~season*reach, random=~1|stream, 
           weights=varComb(vf1,vf2,vf3), data=fall.spring, 
           na.action=na.omit)
anova(M.nds)
summary(M.nds)

#interaction plot
op<-par(mfrow=c(1,1), mar=c(4,4,2,1))
x<-fall.spring[complete.cases(fall.spring),]#strips NA rows from data for plotting
x<-droplevels(x)#removes levels with no values, such as summer
with(x, 
     interaction.plot(reach,season,nrr, 
                      ylim=c(0,8),lty=c(1,12),lwd=2,ylab="NRR", 
                      xlab="Reach", trace.label="Season"))
par(op)
#fall has stronger NRR than spring, 
#but in spring, buried have greater response to added carbon
#whereas in fall, daylight reaches have greater response to added carbon

## Lines 71-79 can be replaced with 81-85
# data = data.frame(aggregate(nrr~season+reach, data=fall.spring, 
#             FUN=function(x) c(mean=mean(x), sd=sd(x), 
#             n=length(x), se=sd(x)/sqrt(length(x)))))
# data#look at aggregated data
# x=data.frame(reach=factor(c("Buried","Buried","Open","Open")), 
#             season=factor(c("Fall","Spring","Fall","Spring"), 
#             levels=c("Fall","Spring")), 
#             nrr.mean=c(4.19558058,2.51384505,7.17557390,1.94443149),
#             nrr.se=c(0.27504630,0.08000855,0.78039976,0.09377607))

x <- group_by(fall.spring, season, reach) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(nrr.mean = mean(nrr, na.rm = TRUE), # na.rm = TRUE to remove missing values
            nrr.sd=sd(nrr, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(nrr)), # of observations, excluding NAs. 
            nrr.se=nrr.sd/sqrt(n))

p<-ggplot(data=x, 
          aes(x=season, y=nrr.mean, fill=reach)) + 
  geom_bar(stat="identity", 
           position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=nrr.mean, 
                    ymax=nrr.mean+nrr.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","snow"))
p+xlab("Season")+ylab("NRR (treatment:control)")+labs(fill="Reach")+theme_bw()
############################################################################################
#Reach scale data analysis


############multivariate approach to incorporate glucose, arabinose, and cellobiose responses
############into a single model
install.packages("vegan")
library(vegan)

#http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html
y <- cbind(reach$glucose.nrr,reach$arabinose.nrr,reach$cellobiose.nrr)
L.y <- cbind(log(reach$glucose.nrr+1),log(reach$arabinose.nrr+1),log(reach$cellobiose.nrr+1))
adonis(y~CBOM.DM+FBOM.DM, data=reach, na.action=na.omit, perm=1e5)
  #runs a simulation that returns CBOM as significant (p=0.039) and FBOM as marginal (p=0.053)
M1<-adonis(y~CBOM.DM+FBOM.DM, data=reach, na.action=na.omit, perm=1e5)
coefficients(M1)#glucose is term 1, arabinose is term 2, cellobiose is term 3
summary(M1)

#interpretation:  the more CBOM and FBOM, the greater the response to added carbon in NDS

plot(reach$glucose.nrr~reach$CBOM.DM)
plot(reach$arabinose.nrr~reach$CBOM.DM)
plot(reach$cellobiose.nrr~reach$CBOM.DM)

#these look terrible...maybe it's better not to graph these?

############N acquisition enzymes
vf1<-varIdent(form = ~1|stream)

M.nace<-gls(NACE.DM~season, weights=vf1, 
         data=reach, na.action=na.omit, method="REML")
anova(M.nace)
summary(M.nace)

model.matrix.gls <- function(M.nace, ...){
  model.matrix(terms(M.nace), data = getData(M.nace), ...)  
}
model.frame.gls <- function(M.nace, ...){
  model.frame(formula(M.nace), data = getData(M.nace), ...)  
}
terms.gls <- function(M.nace, ...){
  terms(model.frame(M.nace),...)  
}

multCompTukey <- glht(M.nace, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#interpretation:  more N acquisition effort in fall compared to spring

data = data.frame(aggregate(NACE.DM~season, data=reach, 
             FUN=function(x) c(mean=mean(x), sd=sd(x), 
             n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(season=factor(c("Fall","Spring","Summer"), 
             levels=c("Fall","Spring","Summer")), 
             NACE.mean=c(1292.0847,473.3486,1020.8635),
             NACE.se=c(154.7755,210.6108,138.3992))
p<-ggplot(data=x, 
          aes(x=season, y=NACE.mean)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=NACE.mean, 
          ymax=NACE.mean+NACE.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black"))

p+xlab("Season")+ylab("NACE (nmol g-1 DM h-1)")+theme_bw()
#check units with Brian


############DOPAH2, correlates to lignin degradation so a metric of recalcitrant carbon use
vf3<-varIdent(form=~1|season)

M.dopah2<-gls(DOPAH2.DM~reach, 
              weights=vf3, data=reach, na.action=na.omit)
anova(M.dopah2)
summary(M.dopah2)

#interpretation:  more lignin degradation effort in buried compared to daylight regardless of season

data = data.frame(aggregate(DOPAH2.DM~reach, data=reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(reach=factor(c("Buried","Daylight"), 
          levels=c("Buried","Daylight")), 
          DOPAH2.mean=c(66.079513,29.633145),
          DOPAH2.se=c(8.583403,9.577104))
p<-ggplot(data=x, 
          aes(x=reach, y=DOPAH2.mean)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=DOPAH2.mean, 
          ymax=DOPAH2.mean+DOPAH2.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black"))

p+xlab("Reach")+ylab("DOPAH2 (nmol g-1 DM h-1)")+theme_bw()
#check units with Brian and write "DOPAH2" properly...wth is this stuff?

############POX, Brian Hill's alternate DOPA metric?
reach$L.POX<-log(reach$POX+1)

M.pox<-gls(L.POX~reach, data=reach, na.action=na.omit)
anova(M.pox)
summary(M.pox)

#interpretation:  more recalcitrant C degradation effort in buried compared to daylight 
    #regardless of season

data = data.frame(aggregate(POX~reach, data=reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(reach=factor(c("Buried","Daylight"), 
          levels=c("Buried","Daylight")), 
          POX.mean=c(753199.1,324964.2),
          POX.se=c(136687.6,163517.4))
p<-ggplot(data=x, 
          aes(x=reach, y=POX.mean)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=POX.mean, 
          ymax=POX.mean+POX.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black"))

p+xlab("Reach")+ylab("POX (mol g-1 C h-1)")+theme_bw()

############LCI, index of carbon lability calculcated as recalcitrant/(labile+recalcitrant)
    #so a smaller number means there is more labile carbon and a bigger number means more
    #recalcitrant carbon

M.lci<-gls(LCI~season+reach+stream, data=reach, na.action=na.omit)
anova(M.lci)
summary(M.lci)

model.matrix.gls <- function(M.lci, ...){
  model.matrix(terms(M.lci), data = getData(M.lci), ...)  
}
model.frame.gls <- function(M.lci, ...){
  model.frame(formula(M.lci), data = getData(M.lci), ...)  
}
terms.gls <- function(M.lci, ...){
  terms(model.frame(M.lci),...)  
}

multCompTukey <- glht(M.lci, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#spring and fall don't differ, but summer has more recalcitrant carbon than the fall
#daylight has more labile carbon than buried

data = data.frame(aggregate(LCI~season+reach, data=reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(reach=factor(c("Buried","Buried","Buried","Open","Open","Open")), 
          season=factor(c("Fall","Spring","Summer","Fall","Spring","Summer"), 
          levels=c("Fall","Spring","Summer")), 
          LCI.mean=c(0.75679360,0.91763295,0.93085557,0.60661323,0.60207482,0.89880824),
          LCI.se=c(0.10750636,0.04472784,0.01133049,0.11575625,0.11305799,0.06472929))
p<-ggplot(data=x, 
          aes(x=season, y=LCI.mean, fill=reach)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=LCI.mean, 
          ymax=LCI.mean+LCI.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black","snow"))
p+xlab("Season")+ylab("LCI (recalcitrant/(labile+recalcitrant))")+labs(fill="Reach")+theme_bw()

############HIX, humification index from Pennino

new.reach<-reach[reach$HIX>0.71,]#subset data to remove outlier

vf1<-varIdent(form = ~1|stream)

M.hix<-gls(HIX~season+reach+stream, 
         weights=vf1, data=new.reach,na.action=na.omit) 
anova(M.hix)
summary(M.hix)

model.matrix.gls <- function(M.hix, ...){
  model.matrix(terms(M.hix), data = getData(M.hix), ...)  
}
model.frame.gls <- function(M.hix, ...){
  model.frame(formula(M.hix), data = getData(M.hix), ...)  
}
terms.gls <- function(M.hix, ...){
  terms(model.frame(M.hix),...)  
}

multCompTukey <- glht(M.hix, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#Fall has higher humification index than spring or summer
#Daylight has higher humification index than buried, possibly due to terrestrial inputs to
  #open reaches

data = data.frame(aggregate(HIX~season+reach, data=new.reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(reach=factor(c("Buried","Buried","Buried","Open","Open","Open")), 
          season=factor(c("Fall","Spring","Summer","Fall","Spring","Summer"), 
          levels=c("Fall","Spring","Summer")), 
          HIX.mean=c(0.889640833,0.875041500,0.894589833,0.938477000,0.881083500,0.874557750),
          HIX.se=c(0.024981359,0.001638500,0.005782358,0.009319571,0.014257931,0.014645250))
p<-ggplot(data=x, 
          aes(x=season, y=HIX.mean, fill=reach)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=HIX.mean, 
          ymax=HIX.mean+HIX.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black","snow"))
p+xlab("Season")+ylab("Humification Index (HIX))")+labs(fill="Reach")+theme_bw()

############BIX, biological freshness index from Pennino

M.bix<-gls(BIX~season, 
              data=reach,na.action=na.omit)  
anova(M.bix)
summary(M.bix)

model.matrix.gls <- function(M.bix, ...){
  model.matrix(terms(M.bix), data = getData(M.bix), ...)  
}
model.frame.gls <- function(M.bix, ...){
  model.frame(formula(M.bix), data = getData(M.bix), ...)  
}
terms.gls <- function(M.bix, ...){
  terms(model.frame(M.bix),...)  
}

multCompTukey <- glht(M.bix, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#spring and summer have higher freshness than fall, but they don't differ from each other

data = data.frame(aggregate(BIX~season, data=reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(season=factor(c("Fall","Spring","Summer"), 
          levels=c("Fall","Spring","Summer")), 
          BIX.mean=c(0.0981945833,0.7412309167,0.7377986000),
          BIX.se=c(0.0007789257,0.0059971767,0.0062660596))
p<-ggplot(data=x, 
          aes(x=season, y=BIX.mean)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=BIX.mean, 
          ymax=BIX.mean+BIX.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black"))

p+xlab("Season")+ylab("Biological Freshness Index (BIX)")+theme_bw()

############FI is fluorescence index from Pennino. bigger=microbial, lower=terrestrial

vf2<-varIdent(form=~1|reach)

M.fi<-gls(FI~season, weights=vf2, 
          data=reach, na.action=na.omit)

anova(M.fi)
summary(M.fi)

model.matrix.gls <- function(M.fi, ...){
  model.matrix(terms(M.fi), data = getData(M.fi), ...)  
}
model.frame.gls <- function(M.fi, ...){
  model.frame(formula(M.fi), data = getData(M.fi), ...)  
}
terms.gls <- function(M.fi, ...){
  terms(model.frame(M.fi),...)  
}

multCompTukey <- glht(M.fi, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#fall has a larger terrestrial signature than spring or summer.  spring and summer not different
  #from each other

data = data.frame(aggregate(FI~season, data=reach, 
          FUN=function(x) c(mean=mean(x), sd=sd(x), 
          n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(season=factor(c("Fall","Spring","Summer"), 
          levels=c("Fall","Spring","Summer")), 
          FI.mean=c(0.335320000,1.278229250,1.262672300),
          FI.se=c(0.028975281,0.004810842,0.008548331))
p<-ggplot(data=x, 
          aes(x=season, y=FI.mean)) + geom_bar(stat="identity", 
          position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=FI.mean, 
          ymax=FI.mean+FI.se), width=0.2, 
          position=position_dodge(0.9)) + scale_fill_manual(values=c("black"))

p+xlab("Season")+ylab("Fluoroscence Index (FI)")+theme_bw()

############P2H is protein to humic ratio from Pennino. bigger numbers are more protein like OM 
  #compared to humic-like OM

vf6<-varPower(form = ~ fitted(.))

M.p2h<-gls(P2H~season+reach+stream, 
        weights=vf6, data=new.reach,na.action=na.omit)

anova(M.p2h)
summary(M.p2h)

model.matrix.gls <- function(M.p2h, ...){
  model.matrix(terms(M.p2h), data = getData(M.p2h), ...)  
}
model.frame.gls <- function(M.p2h, ...){
  model.frame(formula(M.p2h), data = getData(M.p2h), ...)  
}
terms.gls <- function(M.p2h, ...){
  terms(model.frame(M.p2h),...)  
}

multCompTukey <- glht(M.p2h, linfct = mcp(season = "Tukey")) 
summary(multCompTukey)

#spring and summer have higher P2H than fall, daylight has lower P2H than buried with small
  #effect size, possibly due to greater biological activity using the better OM?

data = data.frame(aggregate(P2H~season+reach, data=new.reach, 
        FUN=function(x) c(mean=mean(x), sd=sd(x), 
        n=length(x), se=sd(x)/sqrt(length(x)))))
data#look at aggregated data
x=data.frame(reach=factor(c("Buried","Buried","Buried","Open","Open","Open")), 
        season=factor(c("Fall","Spring","Summer","Fall","Spring","Summer"), 
        levels=c("Fall","Spring","Summer")), 
        P2H.mean=c(0.04535483,0.79201250,0.60917967,0.03616317,0.67722583,0.67767175),
        P2H.se=c(0.01262531,0.06582000,0.04264401,0.01879541,0.08523015,0.03721975))
p<-ggplot(data=x, 
        aes(x=season, y=P2H.mean, fill=reach)) + geom_bar(stat="identity", 
        position=position_dodge(), colour="black") + geom_errorbar(aes(ymin=P2H.mean, 
        ymax=P2H.mean+P2H.se), width=0.2, 
        position=position_dodge(0.9)) + scale_fill_manual(values=c("black","snow"))
p+xlab("Season")+ylab("Protein to Humic Ratio)")+labs(fill="Reach")+theme_bw()
