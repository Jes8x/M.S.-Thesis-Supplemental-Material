library("picante")

#####

##_____________NMDS__##

#Get data in 

data=read.csv(choose.files(), header = TRUE, row.names = 1)
metadata=read.csv(file.choose(), header=TRUE, row.names=1)

all.equal(rownames(metadata),rownames(data))

## Setting seed before running

set.seed(42);ord = metaMDS(data, distance = "euclidean")


#####

##_____________Hclust

dat.dist<-vegdist(data, method = "euclidean") 
dat.clust<- hclust(dat.dist, method = "average")
plot(dat.clust, ylab = "distance")


#####

##_____________Plot__##

dev.new()

#add extra space to the right of the plot
par(mar=c(5, 4, 4, 5), xpd=TRUE)

plot(ord$points, type = "n")
points(ord, choices=c(1,2), display="sites",col=(metadata$Pop),cex=2,pch=as.numeric(metadata$Pop))
legend(cex=1.2, "topright", inset=c(-.15, 0),legend=as.character(paste(" ",unique(metadata$Label))),pch=(unique(metadata$Pop)), col=(unique(metadata$Pop)))


#####
## Plotting vectors

morph <- read.csv(choose.files(), header = TRUE, row.names = 1)

mm <- morph[rownames(morph) %in% rownames(data),]
all.equal(rownames(mm),rownames(data))

colors <- read.csv(choose.files(), header = TRUE, row.names = 1)

fit <- envfit(ord~Avg_Tent_L+Blue, mm)

plot(fit, cex =1)
