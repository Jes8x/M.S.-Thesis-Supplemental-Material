library("vegan")


#### CCA Environment - this is the model with all factors and investigates significance and overlap in variance 

rad_data <- read.csv(file.choose(), header = TRUE, row.names = 1)

morph_data <- read.csv(file.choose(), header = TRUE, row.names = 1)

color_data <- morph_data[rownames(morph_data) %in% rownames(rad_data),]

color_data

set.seed(42); Rad_morph.cca<-cca(rad_data~ Avg_Mouth_W+Avg_Tent_L+Avg_Disc_W+Avg_Num_Tent+Black+Brown+Blue+White+Green+Yellow+Red, data = color_data)
Rad_morph.cca

## This looks at how independent factors overlap
#variance inflation factor VIF >10 is high correlation, more conservative is 2.5 (not good)

vif.cca(Rad_morph.cca)
eigenvals(Rad_morph.cca)


#####
#Checking one variable at a time in this model

# to build our model we next build an empty cca
set.seed(42); lwr.cca <- cca(rad_data ~ 1, data = color_data)
set.seed(42)

## We now run our model starting with the empty one and moving towards the full one.  
#may need to remove parallel option

set.seed(42); modsR2.temp<- ordiR2step(lwr.cca, scope = formula(Rad_morph.cca))
(modsR2.temp)


## Testing the significance of the CCA model:
anova.cca(modsR2.temp)

## Testing the significance of terms (environmental variables):
anova.cca(modsR2.temp, by="terms")

##Testing the significance of CCA axes (at least the first two or three should present a significant p value):
anova.cca(modsR2.temp, by="axis")

# retest for overlap of independent variables
vif.cca(modsR2.temp)

# check for significance of model
temp.permutest<- permutest(modsR2.temp, permutations = how(nperm=999), by = "terms", model = "reduced",parallel = getOption("mc.cores"))
temp.permutest

modsR2.temp$anova
r2adj.temp<-RsquareAdj(modsR2.temp)
r2adj.temp
