library(pavo)
library(colordistance)
library(devtools)
library(recolorize)
library(svMisc)
library(png)


######### Credit to Dr. Jessica H. Arbour for creating the loops
######### within this code and many thanks for her guidance and assistance
#####
#Single image ColorMap

img <- choose.files()

init_fit <- recolorize(img, method = "hist", bins = 2, color_space = "sRGB")


refined_fit <- recluster(init_fit, cutoff = 45)


final_fit <- editLayer(refined_fit, 3, operation = "fill", px_size = 4)



##### Batch stuff?

# get all  images:
images <- dir(pattern = "png", full.names = TRUE)

# make an empty list to store the results:
rc_list <- vector("list", length = length(images))

# run `recolorize2` on each image
# you would probably want to add more sophisticated steps in here as well, but you get the idea

for (i in 1:length(images)) {
  rc_list[[i]] <- suppressMessages(recolorize2(images[i], bins = 2, 
                                               cutoff = 30, plotting = FALSE))
}

# plot for comparison:
layout(matrix(1:4, nrow = 1))
par(mar = rep(0, 4))
for (i in rc_list) {
  plotImageArray(i$original_img)
  plotImageArray(recoloredImage(i))
}

#####
##################################################################
########################color maps################################
##################################################################

#### set to folder with images
images <- dir(pattern = "png", full.names=TRUE)
new.images<-vector(length=length(images),mode="list")

dev.new()

for (i in 1:length(images)) {
  progress(i,max.value=length(images))
  # get an initial fit with generic clustering
  if(i==1){
    
    init_fit <- recolorize(images[i], method = "hist", bins = 3, color_space = "sRGB")
    images# cluster similar colors and fit again
    refined_fit <- recluster(init_fit, cutoff = 50)
    
    
    new.images[[1]]<-refined_fit
  }
  if(i>1){
    new.images[[i]]<-imposeColors(images[i],refined_fit$centers,adjust_centers = FALSE)
  }
  
}

#####

#Exporting

# export color map

for(i in 1:length(new.images)){
recolorize_to_png(new.images[[i]], file=paste0("layer_", i, ".png"))
}

# for (i in 1:length(new.images)) {
#   png::writePNG(new.images[[i]],
#                 target = paste0("layer_", i, ".png"))
# }



#####
#Getting proportion of colors out into csv

library(svMisc)
prop.matrix<-matrix(ncol=7,nrow=16)

for(i in 1:16){
  progress(i,max.value=16)
  props<-table(new.images[[(i+1)]]$pixel_assignments)
  prop.matrix[i,as.numeric(row.names(props[-1]))]<-props[-1]/sum(props[-1])

}

write.csv(prop.matrix, file= "Zoa_color_proportions.csv")


