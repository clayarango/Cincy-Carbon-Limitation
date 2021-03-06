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
install.packages("gridExtra")
library(nlme)
library(nortest)
library(dplyr)
library(multcomp)
library(MASS)
library(ggplot2)
library(gridExtra)

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

reach$season <- factor(reach$season, levels = c("Summer", "Autumn", "Spring"))

############################################################################################
#Carbon limitation analysis with NDS data
fall<-subset(nds, season=="Autumn")
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

x <- group_by(fall.spring, season, reach) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
    summarize(nrr.mean = mean(nrr, na.rm = TRUE), # na.rm = TRUE to remove missing values
         nrr.sd=sd(nrr, na.rm = TRUE),  # na.rm = TRUE to remove missing values
         n = sum(!is.na(nrr)), # of observations, excluding NAs. 
         nrr.se=nrr.sd/sqrt(n))

ggplot(data=x, 
       aes(x=season, y=nrr.mean, fill=reach)) + 
       geom_bar(stat="identity", position=position_dodge(), color = "black") + 
       geom_errorbar(aes(ymin=nrr.mean, ymax=nrr.mean+nrr.se), width=0.2, 
                  position=position_dodge(0.9)) + 
       scale_fill_manual(values=c("black","white")) +
       xlab("Season") +
       ylab("NRR (treatment:control)") +
       ylim(0,8.5) +
       labs(fill="Reach") +
       theme_bw() +
       theme(panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             legend.title=element_text(size=6),
             legend.key=element_blank(),
             legend.position=c(0.5,0.95),
             legend.text=element_text(size=8),
             legend.background=element_blank(),
             legend.direction="horizontal",
             legend.key.size=unit(0.3, "cm"),
             axis.title.y=element_text(size=8),
             axis.title.x=element_text(size=8),
             axis.text.x=element_text(size=8)) +
  annotate("text", x=1.75, y=7, label="season*reach,", size=3, hjust=0) +
  annotate("text", x=1.75, y=6.5, label="p<<0.001", size=3, hjust=0)

ggsave('output/figures/nrrByReachSeason.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

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

plot(reach$glucose.nrr~reach$CBOM.DM, ylab="NRR (treatment:control)", xlab="Coarse Benthic Organic Matter (gAFDM m^2)", pch=19)
plot(reach$arabinose.nrr~reach$CBOM.DM)
plot(reach$cellobiose.nrr~reach$CBOM.DM)
points(reach$arabinose.nrr~reach$CBOM.DM, pch=21)
points(reach$cellobiose.nrr~reach$CBOM.DM, pch=21)
#these look terrible...maybe it's better not to graph these?


nrr.cbom=read.table(file="cincy.nrr.cbom.csv", header=T, sep=",")
nrr.cbom$carbon <- factor(nrr.cbom$carbon, levels = c("Glucose", "Arabinose", "Cellobiose"))
names(nrr.cbom)
label=paste("r^2 == 0.32")
ggplot(nrr.cbom, aes(x=cbom, y=nrr)) +
  geom_errorbar(data=nrr.cbom, mapping=aes(ymin=nrr- nrr.se, ymax=nrr+nrr.se), width=0.2) +
  geom_point(size=2, pch=21, aes(fill=carbon)) +
  scale_fill_manual(name="Carbon", values=c("white", "gray", "black")) +
  xlab(expression(CBOM~gDM~m^{-2})) +
  ylab("NRR (treatment:control)") +
  stat_smooth(method="lm", se=F, col="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.title = element_text(size = 6),
        legend.key = element_blank(),  
        legend.position = c(0.2, 0.85),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(), 
        legend.key.size = unit(0.3, "cm")) +
  annotate("text", x=70, y=18, label="p=0.039", size=3, hjust=0) +
  annotate("text", x=70, y=17, label=label, size=3, parse=T, hjust=0) 

ggsave('output/figures/nrrVcbom.tiff',  # this function is not necessary if creating a 2 panel fig.  see below
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

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

x <- group_by(reach, season) %>%  
  summarize(NACE.mean = mean(NACE.DM, na.rm = TRUE), 
            NACE.sd=sd(NACE.DM, na.rm = TRUE),  
            n = sum(!is.na(NACE.DM)),  
            NACE.se=NACE.sd/sqrt(n))

ggplot(data=x,aes(x=season, y=NACE.mean)) + 
    geom_bar(stat="identity", position=position_dodge(), color = "black") + 
    geom_errorbar(aes(ymin=NACE.mean, ymax=NACE.mean+NACE.se), width=0.2, 
         position=position_dodge(0.9)) + 
    scale_fill_manual(values=c("black")) +
    xlab("Season") +
    ylab(expression(NACE~(nmol~gDM^{-1}~h^{-1}))) + # code for superscripts
    theme_bw() +
    theme(panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         axis.title.y=element_text(size=8),
         axis.title.x=element_text(size=8),
         axis.text.x=element_text(size=8))

ggsave('output/figures/naceBySeason.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#check units with Brian


############DOPAH2, correlates to lignin degradation so a metric of recalcitrant carbon use
vf3<-varIdent(form=~1|season)

M.dopah2<-gls(DOPAH2.DM~reach, 
              weights=vf3, data=reach, na.action=na.omit)
anova(M.dopah2)
summary(M.dopah2)

#interpretation:  more lignin degradation effort in buried compared to daylight regardless of season

x <- group_by(reach, reach) %>%  
  summarize(DOPAH2.mean = mean(DOPAH2.DM, na.rm = TRUE), 
            DOPAH2.sd=sd(DOPAH2.DM, na.rm = TRUE),  
            n = sum(!is.na(DOPAH2.DM)),  
            DOPAH2.se=DOPAH2.sd/sqrt(n))

p.dopa <- ggplot(data=x,aes(x=reach, y=DOPAH2.mean)) + # assign to object to include in two panel fig with POX
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=DOPAH2.mean, ymax=DOPAH2.mean+DOPAH2.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black")) +
  xlab("Reach") +
  ylab(expression(DOPAH2~(nmol~gDM^{-1}~h^{-1}))) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8))

#altered code that doesn't plot x axis label for stacked graph
p.dopa <- ggplot(data=x,aes(x=reach, y=DOPAH2.mean)) + # assign to object to include in two panel fig with POX
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=DOPAH2.mean, ymax=DOPAH2.mean+DOPAH2.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black")) +
  ylab(expression(DOPAH2~(nmol~gDM^{-1}~h^{-1}))) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x=2.3, y=73, label="A", size=6)


# ggsave('output/figures/dopah2ByReach.tiff',  # this function is not necessary if creating a 2 panel fig.  see below
#        units="in",
#        width=3.25,
#        height=3.25,
#        dpi=1200,
#        compression="lzw")

#check units with Brian and write "DOPAH2" properly...wth is this stuff?
#no idea.  I sure hope Brian knows what he is talking about!

############POX, Brian Hill's alternate DOPA metric?
reach$L.POX<-log(reach$POX+1)

M.pox<-gls(L.POX~reach, data=reach, na.action=na.omit)
anova(M.pox)
summary(M.pox)

#interpretation:  more recalcitrant C degradation effort in buried compared to daylight 
    #regardless of season

x <- group_by(reach, reach) %>%  
  summarize(POX.mean = mean(POX, na.rm = TRUE)/1000, # Divide by 1000 to make numbers more manageable
            POX.sd=sd(POX, na.rm = TRUE)/1000,  # Divide by 1000 to make numbers more manageable
            n = sum(!is.na(POX)),  
            POX.se=POX.sd/sqrt(n))

p.pox <- ggplot(data=x,aes(x=reach, y=POX.mean)) + # assign to object to include in two panel fig with POX
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=POX.mean, ymax=POX.mean+POX.se), 
                width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black")) +
  xlab("Reach") +
  ylab(expression(POX~(mmol~gDM^{-1}~h^{-1}))) + #changed to mmol.  -1 should come after DM, not g, right?
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=2.3, y=850, label="B", size=6)

# ggsave('output/figures/poxByReach.tiff', # this function is not necessary if creating a 2 panel fig.  see below
#        units="in",
#        width=3.25,
#        height=3.25,
#        dpi=1200,
#        compression="lzw")

############Two panel plot of DOPA and POX
# Stacked two panel graph.  This make sure left and right edges of plots are aligned.
# We can also do horizontal plot if you prefer.
# Code stolen from http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
gA <- ggplotGrob(p.dopa)  # set up figure
gB <- ggplotGrob(p.pox)  # set up figure
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure

tiff(filename = 'output/figures/dopa.pox.2panel.tiff', #open plotting device
     width = 3.25,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, ncol=1)  # push plot to device
dev.off()  # close device


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

multCompTukey <- glht(M.lci, linfct = mcp(reach = "Tukey")) 
summary(multCompTukey)

#spring and fall don't differ, but summer has more recalcitrant carbon than the fall
#daylight has more labile carbon than buried

x <- group_by(reach, season) %>%  
  summarize(LCI.mean = mean(LCI, na.rm = TRUE), 
            LCI.sd=sd(LCI, na.rm = TRUE),  
            n = sum(!is.na(LCI)),  
            LCI.se=LCI.sd/sqrt(n))

ggplot(data=x, aes(x=season, y=LCI.mean, fill=reach)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_errorbar(aes(ymin=LCI.mean, ymax=LCI.mean+LCI.se), width=0.2, position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +  # consider adopting pure black and white, or a greyscale.  "snow" will likely trigger additiona publication charges
  xlab("Season")+
  ylab("LCI (recalcitrant/(labile+recalcitrant))") +
  ylim(0, 1.05) + # add a bit of room at top for legend
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Eliminate major gridlines
        panel.grid.minor = element_blank(),  # Eliminate minor gridlines
        legend.title = element_text(size = 6),  # Eliminate legend title
        legend.key = element_blank(),  # Eliminate boxes around legend elements
        legend.position = c(0.5, 0.95),  # Specify legend position
        legend.text=element_text(size=8),  # Specify legend text size.  also see legend.title
        legend.background = element_blank(),  # Eliminate white fill in legend.
        legend.direction = "horizontal", # horizontal legend, default is vertical
        legend.key.size = unit(0.3, "cm"), # size of boxes in legend
        axis.title.y = element_text(size = 8), # y axis label text size
        axis.text.y = element_text(size = 8), # y axis tick label text size
        axis.title.x = element_text(size = 8), # x axis label text size
        axis.text.x = element_text(size = 8)) # x axis tick label text size
  
ggsave('output/figures/lciByReachSeason.tiff',  # export as .tif
units="in",  # specify units for dimensions
width=3.25,   # 1 column
height=3.25, # Whatever works
dpi=1200,   # ES&T. 300-600 at PLOS One,
compression = "lzw")

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

x <- group_by(new.reach, season, reach) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(HIX.mean = mean(HIX, na.rm = TRUE), # na.rm = TRUE to remove missing values
            HIX.sd=sd(HIX, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(HIX)), # of observations, excluding NAs. 
            HIX.se=HIX.sd/sqrt(n))

ggplot(data=x, aes(x=season, y=HIX.mean, fill=reach)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") + 
  geom_errorbar(aes(ymin=HIX.mean, ymax=HIX.mean+HIX.se), width=0.2, position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +  
  xlab("Season")+
  ylab("HIX") +
  ylim(0, 1.05) + 
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        legend.title = element_text(size = 6),  
        legend.key = element_blank(),  
        legend.position = c(0.5, 0.95),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(),  
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 8)) 

ggplot(data=new.reach, aes(x=season, y=HIX, fill=reach)) + # consider boxplot for this one
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +  
  xlab("Season")+
  ylab("HIX") +
  ylim(0.84, 0.97) + 
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        legend.title = element_text(size = 6),  
        legend.key = element_blank(),  
        legend.position = c(0.5, 0.95),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(),  
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 8)) 

ggsave('output/figures/box.hixByReachSeason.tiff',  
       units="in",  
       width=3.25,   
       height=3.25, 
       dpi=1200,   
       compression = "lzw")

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

x <- group_by(reach, season) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(BIX.mean = mean(BIX, na.rm = TRUE), # na.rm = TRUE to remove missing values
            BIX.sd=sd(BIX, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(BIX)), # of observations, excluding NAs. 
            BIX.se=BIX.sd/sqrt(n))

x$season <- factor(x$season, levels = c("Summer", "Autumn", "Spring"))

p.BIX <- ggplot(data=x,aes(x=season, y=BIX.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=BIX.mean, ymax=BIX.mean+BIX.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black")) +
  xlab("Season") +
  ylab("BIX") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x=2, y=0.73, label="A", size=6)

ggsave('output/figures/bixBySeason.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

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

x <- group_by(reach, season) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(FI.mean = mean(FI, na.rm = TRUE), # na.rm = TRUE to remove missing values
            FI.sd=sd(FI, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(FI)), # of observations, excluding NAs. 
            FI.se=FI.sd/sqrt(n))

x$season <- factor(x$season, levels = c("Summer", "Autumn", "Spring"))

p.FI <- ggplot(data=x,aes(x=season, y=FI.mean)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=FI.mean, ymax=FI.mean+FI.se), width=0.2, 
                position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black")) +
  xlab("Season") +
  ylab("FI") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8)) +
  annotate("text", x=2, y=1.25, label="B", size=6)

ggsave('output/figures/fiBySeason.tiff',
       units="in",
       width=3.25,
       height=3.25,
       dpi=1200,
       compression="lzw")

#make 2 panel of BIX and FI
gA <- ggplotGrob(p.BIX)  # set up figure
gB <- ggplotGrob(p.FI)  # set up figure
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure

tiff(filename = 'output/figures/bix.fi.2panel.tiff', #open plotting device
     width = 3.25,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, ncol=1)  # push plot to device
dev.off()  # close device




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

multCompTukey <- glht(M.p2h, linfct = mcp(reach = "Tukey")) 
summary(multCompTukey)

#spring and summer have higher P2H than fall, daylight has lower P2H than buried with small
  #effect size, possibly due to greater biological activity using the better OM?

x <- group_by(reach, reach) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(P2H.mean = mean(P2H, na.rm = TRUE), # na.rm = TRUE to remove missing values
            P2H.sd=sd(P2H, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(P2H)), # of observations, excluding NAs. 
            P2H.se=P2H.sd/sqrt(n))

ggplot(data=x, aes(x=season, y=P2H.mean, fill=reach)) + 
  geom_bar(stat="identity", position=position_dodge(), color="black") + 
  geom_errorbar(aes(ymin=P2H.mean, ymax=P2H.mean+P2H.se), width=0.2, position=position_dodge(0.9)) + 
  scale_fill_manual(values=c("black","white")) +  
  xlab("Season")+
  ylab("Protein/Humic ratio") +
  ylim(0, 1.5) + 
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        legend.title = element_text(size = 6),  
        legend.key = element_blank(),  
        legend.position = c(0.5, 0.95),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(),  
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 8)) 

ggsave('output/figures/p2hByReachSeason.tiff',  
       units="in",  
       width=3.25,   
       height=3.25, 
       dpi=1200,   
       compression = "lzw")

########################################################################
#New four panel figure with HIX, BIX, P/H, FI   20-Jun-17
########################################################################

reach<-read.table(file="cincy.reach.csv", header=T, sep=',')

reach$season <- factor(reach$season, levels = c("Spring", "Summer", "Autumn"))

new.reach<-reach[reach$HIX>0.71,]#subset data to remove outlier

install.packages("ggplot2")
install.packages("gridExtra")
install.packages("scales")
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

p.HIX <- ggplot(data=new.reach, aes(x=season, y=HIX, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +  
  xlab("Season")+
  ylab("HIX") +
  ylim(0.84, 0.97) + 
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        legend.title = element_text(size = 6),  
        legend.key = element_blank(),  
        legend.position = c(0.5, 0.95),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(),  
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x=0.8, y=0.96, label="A", size=6) +
  annotate("text", x=1, y=0.91, label="a", size=3) +
  annotate("text", x=2, y=0.91, label="a", size=3) +
  annotate("text", x=3, y=0.955, label="b", size=3) +
  annotate("text", x=1, y=0.94, label="reach,", size=3, hjust=0) +
  annotate("text", x=1, y=0.935, label="p=0.021", size=3, hjust=0)

p.P2H <- ggplot(data=reach, aes(x=season, y=P2H, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +  
  xlab("Season")+
  ylab("Protein/Humic ratio") +
  ylim(0, 1.0) + 
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 8),
        legend.position = "none") +
  annotate("text", x=3, y=0.95, label="C", size=6) +
  annotate("text", x=1, y=0.9, label="a", size=3) +
  annotate("text", x=2, y=0.8, label="a", size=3) +
  annotate("text", x=3, y=0.15, label="b", size=3) +
  annotate("text", x=1, y=0.3, label="reach,", size=3, hjust=0) +
  annotate("text", x=1, y=0.25, label="p=0.0002", size=3, hjust=0)

p.BIX <- ggplot(data=reach,aes(x=season, y=BIX, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +
  xlab("Season") +
  ylab("BIX") +
  ylim(0, 0.8) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position ="none") +
  annotate("text", x=3, y=0.73, label="B", size=6) +
  annotate("text", x=1, y=0.8, label="a", size=3) +
  annotate("text", x=2, y=0.8, label="a", size=3) +
  annotate("text", x=3, y=0.15, label="b", size=3)

p.FI <- ggplot(data=reach,aes(x=season, y=FI, fill=reach)) + 
  #scale_y_log10(breaks = c(0.5, 1), # specify where to have ticks
  #             labels = c(0.5,1)) + # specify tick labels
  geom_boxplot() +
  scale_fill_manual(values=c("black", "white")) +
  xlab("Season") +
  ylab("FI") +
  ylim(0, 1.5) + # this was overriding the log scale
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.position = "none") +
  annotate("text", x=3, y=1.4, label="D", size=6) +
  annotate("text", x=1, y=1.4, label="a", size=3) +
  annotate("text", x=2, y=1.4, label="a", size=3) +
  annotate("text", x=3, y=0.55, label="b", size=3)

#make 4 panel 
gA <- ggplotGrob(p.HIX)  # set up figure
gB <- ggplotGrob(p.BIX)  # set up figure
gC <- ggplotGrob(p.P2H)  # set up figure
gD <- ggplotGrob(p.FI)  # set up figure
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5], gD$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure
gC$widths[2:5] <- as.list(maxWidth)  # set up figure
gD$widths[2:5] <- as.list(maxWidth)  # set up figure

tiff(filename = 'output/figures/hix.bix.p2h.fi.4panel.tiff', #open plotting device
     width = 6.5,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, gC, gD, ncol=2)  # push plot to device
dev.off()  # close device

########################################################################
#New two panel figure with POX and NACE (dropping DOPAH2)   20-Jun-17
########################################################################

p.POX <- ggplot(data=reach, aes(x=season, y=POX/1000, fill=reach)) + # assign to object to include in two panel fig with POX
  geom_boxplot() + 
  scale_fill_manual(values=c("black", "white")) +
  xlab("Season") +
  ylab(expression(POX~(mmol~gDM^{-1}~h^{-1}))) + #changed to mmol.  -1 should come after DM, not g, right?
  ylim(0, 1600) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        legend.title = element_text(size = 6),  
        legend.key = element_blank(),  
        legend.position = c(0.5, 0.95),  
        legend.text=element_text(size=8),  
        legend.background = element_blank(),  
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3, "cm"), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x=3, y=1500, label="A", size=6) +
  annotate("text", x=2.65, y=1300, label="reach,", size=3, hjust=0) +
  annotate("text", x=2.65, y=1200, label="p=0.004", size=3, hjust=0)

p.NACE <- ggplot(data=reach, aes(x=season, y=NACE.DM, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black", "white")) +
  xlab("Season") +
  ylab(expression(NACE~(nmol~gDM^{-1}~h^{-1}))) + # code for superscripts
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.position ="none") +
  annotate("text", x=0.8, y=1900, label="B", size=6) +
  annotate("text", x=1, y=1200, label="a", size=3) +
  annotate("text", x=2, y=1600, label="b", size=3) +
  annotate("text", x=3, y=1800, label="c", size=3)

############Two panel plot of DOPA and POX
# Stacked two panel graph.  This make sure left and right edges of plots are aligned.
# We can also do horizontal plot if you prefer.
# Code stolen from http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
gA <- ggplotGrob(p.POX)  # set up figure
gB <- ggplotGrob(p.NACE)  # set up figure
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])  # set up figure
gA$widths[2:5] <- as.list(maxWidth)  # set up figure
gB$widths[2:5] <- as.list(maxWidth)  # set up figure

tiff(filename = 'output/figures/pox.nace.2panel.tiff', #open plotting device
     width = 3.25,
     height = 6.5,
     units = "in",
     res = 1200,
     compression = "lzw")
grid.arrange(gA, gB, ncol=1)  # push plot to device
dev.off()  # close device

########################################################################
#New LCI figure...called LCI in the data table, but heading is CQI
#CQI in data table is not correct     20-Jun-17
########################################################################

ggplot(data=reach, aes(x=season, y=LCI, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +  # consider adopting pure black and white, or a greyscale.  "snow" will likely trigger additiona publication charges
  xlab("Season")+
  ylab("CQI (lnPOX/(lnBG+lnPOX))") +
  ylim(0, 1.05) + # add a bit of room at top for legend
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Eliminate major gridlines
        panel.grid.minor = element_blank(),  # Eliminate minor gridlines
        legend.title = element_text(size = 6),  # Eliminate legend title
        legend.key = element_blank(),  # Eliminate boxes around legend elements
        legend.position = c(0.5, 0.10),  # Specify legend position
        legend.text=element_text(size=8),  # Specify legend text size.  also see legend.title
        legend.background = element_blank(),  # Eliminate white fill in legend.
        legend.direction = "horizontal", # horizontal legend, default is vertical
        legend.key.size = unit(0.3, "cm"), # size of boxes in legend
        axis.title.y = element_text(size = 8), # y axis label text size
        axis.text.y = element_text(size = 8), # y axis tick label text size
        axis.title.x = element_text(size = 8), # x axis label text size
        axis.text.x = element_text(size = 8)) +
  annotate("text", x=1, y=1.025, label="ab", size=3) +
  annotate("text", x=2, y=1.025, label="b", size=3) +
  annotate("text", x=3, y=1.025, label="a", size=3) +
  annotate("text", x=2.5, y=0.35, label="reach,", size=3, hjust=0) +
  annotate("text", x=2.5, y=0.25, label="p=0.001", size=3, hjust=0)# x axis tick label text size

ggsave('output/figures/cqiByReachSeason.tiff',  # export as .tif
       units="in",  # specify units for dimensions
       width=3.25,   # 1 column
       height=3.25, # Whatever works
       dpi=1200,   # ES&T. 300-600 at PLOS One,
       compression = "lzw")
