bge = read.table(file="bge.csv", header=T, sep=",")
names(bge)

M1=gls(vf~season+reach+stream, data=bge, na.action=na.omit, method="ML")
summary(M1)
anova(M1)

plot(M1)

M2<-lme(vf~season+reach, random=~1|stream, 
        data=bge, na.action=na.omit, method="REML") 

plot(M2)

anova(M1,M2)

vf1 = varIdent(form = ~ 1|reach)
vf2 = varIdent(form= ~ 1|season)
vf3 = varIdent(form= ~ 1|stream)

M3=lme(vf~season+reach, random=~1|stream, data=bge, weights= vf1, na.action=na.omit, method="REML")
M4=lme(vf~season+reach, random=~1|stream, data=bge, weights= vf2, na.action=na.omit, method="REML")
M5=lme(vf~season+reach, random=~1|stream, data=bge, weights= vf3, na.action=na.omit, method="REML")

anova(M2,M4)

plot(M4)
anova(M4)

x <- group_by(bge, season, reach) %>%  # Grouping function causes subsequent functions to aggregate by season and reach
  summarize(vf.avg = mean(vf, na.rm = TRUE), # na.rm = TRUE to remove missing values
            vf.sd=sd(vf, na.rm = TRUE),  # na.rm = TRUE to remove missing values
            n = sum(!is.na(vf)), # of observations, excluding NAs. 
            vf.se=vf.sd/sqrt(n))

ggplot(data=bge, aes(x=season, y=vf, fill=reach)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("black","white")) +  # consider adopting pure black and white, or a greyscale.  "snow" will likely trigger additiona publication charges
  xlab("Season")+
  ylab(expression(V[f]~~(mm~s^{-1}))) +
  labs(fill="Reach")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Eliminate major gridlines
        panel.grid.minor = element_blank(),  # Eliminate minor gridlines
        legend.title = element_text(size = 6),  # Eliminate legend title
        legend.key = element_blank(),  # Eliminate boxes around legend elements
        legend.position = c(0.5, 0.9),  # Specify legend position
        legend.text=element_text(size=8),  # Specify legend text size.  also see legend.title
        legend.background = element_blank(),  # Eliminate white fill in legend.
        legend.direction = "horizontal", # horizontal legend, default is vertical
        legend.key.size = unit(0.3, "cm"), # size of boxes in legend
        axis.title.y = element_text(size = 8), # y axis label text size
        axis.text.y = element_text(size = 8), # y axis tick label text size
        axis.title.x = element_text(size = 8), # x axis label text size
        axis.text.x = element_text(size = 8)) +
  annotate("text", x=1, y=0.0125, label="reach,", size=3, hjust=0) +
  annotate("text", x=1, y=0.0115, label="p=0.056", size=3, hjust=0)

ggsave('output/figures/vfByReachSeason.tiff',  # export as .tif
       units="in",  # specify units for dimensions
       width=3.25,   # 1 column
       height=3.25, # Whatever works
       dpi=1200,   # ES&T. 300-600 at PLOS One,
       compression = "lzw")
