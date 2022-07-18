library("picante")
library("adegenet")
library("devtools")
library("ggplot2")

#####

####Importing data - Must change file extension manually

obj1 <- import2genind("populations.snps.gen")


#Turning genpop file into csv

sum(is.na(obj1$tab))
Symbio_Genpop <- tab(obj1, freq = TRUE, NA.method = "mean")
write.csv(Symbio_Genpop,file="Symbio_Genpop.csv")


#####

#Mantel test for Sanger and 2b-RAD

#Read tree into object
phylo <-read.tree(choose.files())

phylo

#Read 2b-RAD data into object

comm <- read.csv(choose.files(), header = TRUE, row.names = 1)

comm[1:5, 1:5]

#attributes(comm)

#Use this to automatically reorder rows using tree and community data (i.e. 2b-RAD genpop data)
#check for mismatches/missing species
#using .data instead of .comm since community data is in data.frame class
combined <- match.phylo.data(phylo, comm)

# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy <- combined$phy
comm <- combined$data #this needs to be $data and not $comm since we used the ".data" ending

#attributes(phy)

#write comm data to new csv

write.csv(comm, file = "Zoa_2bRAD_matched_to_tree.csv")


#Get distances between terminal nodes)

Zoa_Tree_Distance <-cophenetic.phylo(phy)

write.csv(Zoa_Tree_Distance, file = "Zoa_Tree_matched_to_RAD.csv")

#Read in our csv files

RawRAD=read.csv(choose.files(),header=TRUE, row.names=1)
Sanger=read.csv(choose.files(),header=TRUE, row.names=1)

#2b-RAD data needs to be transformed into Distance matrix

RAD <- vegdist(RawRAD, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

Sangerveg <- vegdist(Sanger, method = "euclidean", binary = FALSE, diag = TRUE, na.rm = FALSE)

#Note that in the RAD data since it is a "dist" class the sample names are known as "labels" not rownames
all.equal(labels(RAD),labels(Sangerveg))



mantel(RAD, Sanger)




#####

#Comparing Zoa 2b-RAD to Symbio 2b-RAD

Zoa_RAD_Excel <- read.csv(choose.files(), header = TRUE, row.names=1)

Symbio_RAD_Excel <- read.csv(choose.files(), header = TRUE, row.names = 1)

#Trying to remove samples from Zoa data that are not in Symbio data
#Doing this manually since sample names do not match


Zoa_RAD_Excel[1:8, 1:8]

Symbio_RAD_Excel[1:8, 1:8]

Zoa_rm_rows <- Zoa_RAD_Excel[-c(7,8,10,13),]

Zoa_rm_rows[1:16,1:16]

Symbio_veg <- vegdist(Symbio_RAD_Excel, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

Zoa_veg <- vegdist(Zoa_rm_rows, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

Zoa_veg

Symbio_veg

#this wont work since names do not match up
#all.equal(labels(Zoa_veg),labels(Symbio_veg))

#Get both data sets into distance matrix

mantel(Zoa_veg, Symbio_veg)

#####

#Getting morphological data in 

morpho <- read.csv(choose.files(), header = TRUE, row.names=1)

morpho

#Have to match up sample names with tree file


#Mantel Test with Sanger and morpho data

#Read tree into object
phylo_m <-read.tree(choose.files())

phylo_m


#check for mismatches/missing species
#using .data instead of .comm since community data is in data.frame class
combined_m <- match.phylo.data(phylo_m, morpho)

# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy_morpho <- combined_m$phy
morpho_matched_sanger <- combined_m$data #this needs to be $data and not $comm since we used the ".data" ending

#Writing morpho metrics to csv in new order
write.csv(morpho_matched_sanger, file = "morpho_color_matched_to_sanger.csv")

#Writing tree to new csv because it's now matched to morpho data

Zoa_tree_to_morpho <-cophenetic.phylo(phy_morpho)

write.csv(Zoa_tree_to_morpho, file = "tree_matched_to_morpho_color.csv")

#Read in CSVs and check sample order 

morpho_reordered=read.csv(choose.files(),header=TRUE, row.names=1)
Sanger_morpho=read.csv(choose.files(),header=TRUE, row.names=1)

all.equal(row.names(morpho_reordered),rownames(Sanger_morpho))

morpho_dist <- vegdist(morpho_reordered, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

morpho_dist

Sanger_veg <- vegdist(Sanger_morpho, method = "euclidean", binary = FALSE, diag = TRUE, na.rm = FALSE)

Sanger_veg

morph_color_sanger_mantel <- mantel(morpho_dist, Sanger_veg)

morph_color_sanger_mantel



#####

#Mantel test for Morpho and 2b-RAD data

morpho_r <- read.csv(choose.files(), header = TRUE, row.names=1)

RAD_r <- read.csv(choose.files(), header = TRUE, row.names = 1)

#####
#This section is for matching Symbiodinium RAD data to Zoa Morpho

morpho_r

RAD_r[1:16,1:16]

morpho_rm <- morpho_r[-c(4,8,9,11,12,15),]

morpho_rm

morpho_rm_veg <- vegdist(morpho_rm, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

Symbio_RAD_veg <- vegdist(RAD_r, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)


mantel(morpho_rm_veg, Symbio_RAD_veg)


#####
#Section for morpho to Zoa RAD data

#Trying to remove samples from morpho data that are not in RAD data

matched_morpho_r <- morpho_r[rownames(morpho_r) %in% rownames(RAD_r),]

matched_morpho_r

all.equal(rownames(RAD_r),rownames(matched_morpho_r))

#Get both data sets into distance matrix

morpho_r_dist <- vegdist(matched_morpho_r, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

morpho_r_dist

RAD_r_dist <- vegdist(RAD_r, method="euclidean", binary=FALSE, diag = TRUE, na.rm = FALSE)

RAD_r_dist

mantel(morpho_r_dist, RAD_r_dist)


