# SCRIPT TO CREATE REDUCED VERSION OF FIG 3 FROM BEAULIEU ET AL 2014.
# ELIMINATING WINTER SEASON FROM FIG.
# DATA CONTAINED IN aggregated.om.for.Sigma.Plot.txt WHICH WAS CREATED
# DURING THE PREPARATION OF THE 2014 PAPER.  FILE PASTED INTO P ROJECT
# AND READ IN HERE.

# LIBRARIES
library(grid)
library(ggplot2)
library(dplyr)

# Data
df <- read.table("aggregated.om.for.Sigma.Plot.txt", 
                 as.is = TRUE, header = TRUE)
str(df)             

# Strip out winter season
df <- filter(df, season != "WINTER")

legend <- data.frame (x=c(2.3,2.3), y=c(50,45), shape=c(21,22), label=c("OPEN", "Buried"))
PanelA <-
  ggplot(df, aes(season, CBOM)) +   geom_line(aes(group = reach)) + 
  geom_errorbar(aes(ymax = CBOM + CBOM.SE, ymin = CBOM - CBOM.SE), width = 0.1) + 
  scale_shape_identity() +
  geom_point(aes(shape = shape), fill = 'white', size = 3) +
  geom_point(data=legend, aes(x=x, y=y, shape = shape), fill = 'white', size = 3) +
  annotate(geom='text', x=c(0.7, 3, 3), y=c(55, 50, 45), label=c("a","OPEN", "BURIED"), size = 3) +
  theme_bw() +
  ylab(expression(paste('CBOM (g AFDM ', m^-2, ')'))) +
  theme(axis.text.x = element_text(size =6),
        axis.text.y = element_text(size =8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size =8, vjust=0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


PanelB <-
  ggplot(df, aes(season, FBOM)) + geom_line(aes(group = reach)) + 
  geom_errorbar(aes(ymax = FBOM + FBOM.SE, ymin = FBOM - FBOM.SE), width = 0.1) + 
  scale_shape_identity() +
  geom_point(aes(shape = shape), fill = 'white', size = 3) +
  annotate(geom='text', x=0.7, y=30, label="b", size = 3) +
  theme_bw() +
  ylab(expression(paste('FBOM (g AFDM ', m^-2, ')'))) +
  theme(axis.text.x = element_text(size =6),
        axis.text.y = element_text(size =8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size =8, vjust=0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

PanelC <-
  ggplot(df, aes(season, chl)) + geom_line(aes(group = reach)) + 
  geom_errorbar(aes(ymax = chl + chl.se, ymin = chl - chl.se), width = 0.1) + 
  scale_shape_identity() +
  geom_point(aes(shape = shape), fill = 'white', size = 3) +
  annotate(geom='text', x=0.7, y=215, label="c", size = 3) +
  theme_bw() +
  ylab(expression(paste('Chlorophyll a (mg ', m^-2, ')'))) +
  theme(axis.text.x = element_text(size =6),
        axis.text.y = element_text(size =8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size =8, vjust=0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

PanelD <-
  ggplot(df, aes(season, PERI)) + geom_line(aes(group = reach)) + 
  geom_errorbar(aes(ymax = PERI + PERI.SE, ymin = PERI - PERI.SE), width = 0.1) + 
  scale_shape_identity() +
  geom_point(aes(shape = shape), fill = 'white', size = 3) +
  annotate(geom='text', x=0.7, y=35, label="d", size = 3) +
  theme_bw() +
  ylab(expression(paste('Periphyton (g AFDM ', m^-2, ')'))) +
  theme(axis.text.x = element_text(size =6),
        axis.text.y = element_text(size =8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size =8, vjust=0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )   


tiff(filename="omStandingStock.tiff",
     width=4, height=4, units="in",
     res=800,compression="lzw")

grid.newpage()          
pushViewport(viewport(layout = grid.layout(2,2))) 
print(PanelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(PanelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(PanelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(PanelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
