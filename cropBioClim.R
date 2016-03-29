#Crops all the biolcim to north america
library(raster)

file_list=list.files('~/data/bioclim/original/current/', pattern='*bil', full.names = TRUE)

file_list=file_list[!grepl('xml', file_list)]

file_list=sort(file_list)

north_america_extent=extent(-130, -60, 25, 50)

outputDir='./bioclim_cropped/'
for(this_file in file_list){
  rasterObj=raster(this_file)
  this_file_base=basename(this_file)
  rasterObjClipped=crop(rasterObj, north_america_extent, filename=paste(outputDir, this_file_base,sep=''))
}