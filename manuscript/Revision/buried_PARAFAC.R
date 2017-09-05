# Buried PARAFAC

# Working with Data on J:\UMBC\Buried Stream\Revisions

DF.COMP4 = read.csv("J:/UMBC/Buried Stream/Revisions/4Comp.csv")  
DF.COMP3_original = read.csv("J:/UMBC/Buried Stream/Revisions/3Comp.csv")  

# Combine DF.COMP4 DF.COMP3 (just SQ1 and SQ2 componsents)
# Combining 
DF.COMP4_sub = DF.COMP4[,c("ID","Site","Season","Reach","Comp2","Comp4")]
DF.COMP3_sub = DF.COMP3_original[,c("ID","Site","Season","Reach","Comp2","Comp3")]

names(DF.COMP4_sub)[which(colnames(DF.COMP4_sub)=="Comp2")] = "SQ1"
names(DF.COMP4_sub)[which(colnames(DF.COMP4_sub)=="Comp4")] = "SQ2"
names(DF.COMP3_sub)[which(colnames(DF.COMP3_sub)=="Comp2")] = "SQ1"
names(DF.COMP3_sub)[which(colnames(DF.COMP3_sub)=="Comp3")] = "SQ2"

comp_all = rbind(DF.COMP4_sub,DF.COMP3_sub)

##############################################################################
boxplot(SQ1~Reach, data=comp_all,ylab="Humic-terrestiral POM")
boxplot(SQ1~Site, data=comp_all,ylab="Humic-terrestiral POM")
boxplot(SQ1~Season, data=comp_all,ylab="Humic-terrestiral POM")

boxplot(SQ2~Reach, data=comp_all,ylab="Humic DOM")
boxplot(SQ2~Site, data=comp_all,ylab="Humic DOM")
boxplot(SQ2~Season, data=comp_all,ylab="Humic DOM")



#************
# DF.COMP4
boxplot(Comp1~Reach, data=DF.COMP4,ylab="Authochthonous DOM")
boxplot(Comp1~Site, data=DF.COMP4,ylab="Authochthonous DOM")
boxplot(Comp1~Season, data=DF.COMP4,ylab="Authochthonous DOM")

boxplot(Comp2~Reach, data=DF.COMP4,ylab="Humic-terrestiral POM")
boxplot(Comp2~Site, data=DF.COMP4,ylab="Humic-terrestiral POM")
boxplot(Comp2~Season, data=DF.COMP4,ylab="Humic-terrestiral POM")

boxplot(Comp3~Reach, data=DF.COMP4,ylab="Authochthonous DOM")
boxplot(Comp3~Site, data=DF.COMP4,ylab="Authochthonous DOM")
boxplot(Comp3~Season, data=DF.COMP4,ylab="Authochthonous DOM")

boxplot(Comp4~Reach, data=DF.COMP4,ylab="Humic DOM")
boxplot(Comp4~Site, data=DF.COMP4,ylab="Humic DOM")
boxplot(Comp4~Season, data=DF.COMP4,ylab="Humic DOM")

#************
# DF.COMP3
# Remove the one spring site EASBUSPS1

DF.COMP3 <- DF.COMP3[!DF.COMP3$ID == "EASBUSPS1 ",]

boxplot(Comp1~Reach, data=DF.COMP3,ylab="Authochthonous DOM")

boxplot(Comp2~Reach, data=DF.COMP3,ylab="Humic-terrestiral POM")

boxplot(Comp3~Reach, data=DF.COMP3,ylab="Authochthonous DOM")


#############################################################################
t.test(Comp1~Reach, data=DF.COMP4)


# non-parametrics tests of group differences
# All seasons
wilcox.test(comp_all$SQ1~comp_all$Reach) # where y and x are numeric 
# p-value = 0.159

wilcox.test(comp_all$SQ2~comp_all$Reach) # where y and x are numeric 
# p-value = 0.1634

#************************************
# for summer, winter spring
wilcox.test(DF.COMP4$Comp1~DF.COMP4$Reach) # where y and x are numeric 
# p-value = 0.1633

wilcox.test(DF.COMP4$Comp2~DF.COMP4$Reach) # where y and x are numeric 
# p-value = 0.05084

wilcox.test(DF.COMP4$Comp3~DF.COMP4$Reach) # where y and x are numeric 
# p-value = 0.04186

wilcox.test(DF.COMP4$Comp4~DF.COMP4$Reach) # where y and x are numeric 
# p-value = 0.1064

#************************************
# for just Fall
wilcox.test(DF.COMP3$Comp1~DF.COMP3$Reach) # where y and x are numeric 
# p-value = 0.3777

wilcox.test(DF.COMP3$Comp2~DF.COMP3$Reach) # where y and x are numeric 
# p-value = 0.7987

wilcox.test(DF.COMP3$Comp3~DF.COMP3$Reach) # where y and x are numeric 
# p-value = 0.4776



# non-parametric ANOVA

# Seasonal Differences
pairwise.wilcox.test(comp_all$SQ1, comp_all$Season, p.adj = "bonf", exact=F, paired=F)
# For Humic-Terrestrial Winder is diff than summer and spring
pairwise.wilcox.test(comp_all$SQ2, comp_all$Season, p.adj = "bonf", exact=F, paired=F)
# For Humic-Terrestrial Winder is diff than summer and spring
