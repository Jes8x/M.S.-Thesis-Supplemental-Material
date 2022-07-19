library("picante")

phylo <- read.tree(file.choose())
plot(phylo)
morpho <- read.csv(file.choose(), header = TRUE, row.names=1)
rownames(morpho)


#check for mismatches/missing species
#using .data since community data is in data.frame class
combined <- match.phylo.data(phylo, morpho)

# the resulting object is a list with $phy and $comm elements.  replace our
# original data with the sorted/matched data
phy_matched <- combined$phy
x <- combined$data #this needs to be $data and not $comm since we used the ".data" ending

phy_matched

x

#phylo_sig <- phylosignal(x$Avg_Mouth_W, phy_matched, reps = 999)
multiPhylosignal(x, multi2di(phy_matched))

## plot

plot(phy_matched, , show.tip.label = T, cex = 0.7, label.offset = .25)
tiplabels(pch = 19, col = "black", cex = .5 * (x$Avg_Num_Tent))

dev.new()
plot(phy_matched, show.tip.label = TRUE, show.node.label = TRUE, cex = 0.7,label.offset = 2.25)
tiplabels(pch = 19, col = "#000000",cex = 3 * (x[, "Avg_Tent_L"]/max(x[,"Avg_Tent_L"])))
#tiplabels(pch = 19, col = "#009E73",cex = 3 * ((x[, "Avg_Mouth_W"]+22)/(max(x[,"Avg_Mouth_W"]+22))), adj = 1)
#tiplabels(pch = 19, col = "#e79f00",cex = 3 * (x[, "Avg_Disc_W"]/max(x[,"Avg_Disc_W"])), adj = 1.5)
#tiplabels(pch = 19, col = "#9ad0f3",cex = 3 * (x$Avg_Num_Tent/max(x$Avg_Num_Tent)), adj = 2)
tiplabels(pch = 19, col = "#0072B2",cex = 3 * (x$Blue/max(x$Blue)), adj = 2.5) ## set at .54 for fig w/4 metrics
#tiplabels(pch = 19, col = "#D55E00",cex = 2*(x$Blue /max(x$Blue)), adj = 3) ## set at .56 for fig w/4 metrics
#tiplabels(pch = 19, col = "gray",cex = 2 * (x[, "X.N"]/max(x[,"X.N"])), adj = .56)
#tiplabels(pch = 19, col = "#CC79A7",cex = 3 * (x[, "X.C"]/max(x[,"X.C"])), adj = .58)
#tiplabels(pch = 19, col = "yellow3",cex = 3 * (x[, "C.N"]/max(x[,"C.N"])), adj = .60)
