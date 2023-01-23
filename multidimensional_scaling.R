library(ggdendro)

rownames(pairwise_language_distances) = geodata$glotto_code

d = pairwise_language_distances
d = dist(pairwise_language_distances, method = "manhattan")

hc <- d %>%
  as.dist() %>%
  hclust()

affine_tree = as.phylo(hc)

statuses = QuartetStatus(affine_tree, my_tree)
SimilarityMetrics(statuses, similarity=T)

ggdendrogram(hc, rotate = T, size=.05)

nClasses = 28

soundClasses = cutree(hc, k=nClasses) %>% as.numeric()
l = geodata %>% .[[1]]


for (i in 1:nClasses) {
  print(paste0(l[soundClasses == i], collapse=" "))
}


setwd("C:/Users/stein/Dropbox/Studium/7. Semester/BA-Thesis/BA-Thesis_NorthEuraLex")

pairwise_language_distances = read.csv("language_distances.csv")
pairwise_language_distances = read.csv("language_distances_averaged.csv")
pairwise_language_distances = read.csv("ldn.csv")
pairwise_language_distances = read.csv("ldnd.csv")
pairwise_language_distances = read.csv("ldnd_notscaledto1.csv")

library(spData)
library(sf)
library(tidyverse)

geodata = read_tsv("northeuralex-0.9-language-data.tsv", show_col_types = FALSE)
head(geodata)

# apply multi-dimensional scaling on the resulting matrix
pairwise_language_distances = pairwise_language_distances - min(pairwise_language_distances)
pairwise_language_distances = max(pairwise_language_distances) - pairwise_language_distances
# scaled_pairwise_distances = scaled_pairwise_distances / max(scaled_pairwise_distances)
f = as.matrix(pairwise_language_distances)
diag(f) = 0
f_log = f + 1
f_log = log(f_log)


multi_dimensional_scaling_raw = cmdscale(pairwise_language_distances, k=3)

multi_dimensional_scaling = as.tibble(multi_dimensional_scaling_raw)

iso_codes = geodata %>% select(iso_code)

mds = 
  multi_dimensional_scaling %>%
  mutate(iso_codes = iso_codes$iso_code) %>%
  # subtract the smallest value
  mutate(x = V1 - min(.$V1), y = V2 - min(.$V2), z = V3 - min(.$V3)) %>%
  # divide by maximum value
  mutate(r = x / max(.$x), g = y / max(.$y), b = z / max(.$z)) %>%
  # combine r, g, and b into a color
  mutate(col = rgb(r,g,b)) %>%
  # select the only two relevant columns
  select(iso_codes, col, r, g, b) 

head(mds)


# join coordinates with the color data 
colored_geodata = geodata %>% 
  rename(iso_codes = iso_code) %>%
  select(iso_codes, latitude, longitude, family, subfamily) %>%
  inner_join(mds)


# produce an extremely slim polygon that represents the 60° W meridian
meridian_55_west = 
  st_polygon(x = list(rbind(c(-55.0001, 90), # upper left corner of the polygon
                            c(-55, 90), # upper right corner of the polygon
                            c(-55, -90), # lower right corner of the polygon
                            c(-55.0001, -90), # lower left corner of the polygon
                            c(-55.0001, 90)))) %>%
  st_sfc() %>%
  # set the crs of the polygon to the geodetic system for world
  st_set_crs(4326)

# remove everything on the 60° meridian from our world data, effectively cutting all polygons that cross the meridian
world_without_55 = 
  world %>% 
  st_difference(meridian_55_west) %>% 
  st_transform("+proj=eqearth lon_0=125")


# convert the coordinates in geodata into simple features
geodata_sf = 
  colored_geodata %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>%
  st_transform("+proj=eqearth lon_0=125")

library(tmap)

#tmap_save(
geodata_sf %>%
  tm_shape() +
  tm_symbols(size=0.2) +
  tm_shape(world_without_55) +
  tm_polygons(alpha = 1, col = 'white', border.alpha = 0.4) +
  tm_layout(bg.color = 'lightblue') +
  tm_shape(geodata_sf) +
  tm_symbols(size=0.5, col="col", border.lwd=1, alpha=1)#, filename="mapminus1.png" ,dpi = 600)


glotto_codes = as.list(geodata$glotto_code) %>% 
  unlist()
colnames(pairwise_language_distances) = glotto_codes
row.names(pairwise_language_distances) = glotto_codes
head(pairwise_language_distances)









geodata_uralic = geodata_sf %>% filter(subfamily == "Balto-Slavic")

bbox_new = st_bbox(geodata_uralic)
xrange = bbox_new$xmax - bbox_new$xmin
yrange = bbox_new$ymax - bbox_new$ymin

bbox_new[1] <- bbox_new[1] - (0.1 * xrange) # xmin - left
# bbox_new[3] <- bbox_new[3] + (0.1 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.05 * yrange) # ymin - bottom
# bbox_new[4] <- bbox_new[4] + (0.1 * yrange) # ymax - top

bbox_new = bbox_new %>%
  st_as_sfc()

geodata_uralic %>%
  tm_shape(bbox = bbox_new) +
  tm_symbols(size=0.2) +
  tm_shape(world_without_55) +
  tm_polygons(alpha = 1, col = 'white', border.alpha = 0.3) +
  tm_layout(bg.color = 'lightblue') +
  tm_shape(geodata_uralic) +
  tm_symbols(size=0.5, col="col", border.lwd=1, alpha=1) +
  tm_text('subfamily', size = 1/2, just = 'right')



library(glottoTrees)
my_list = list(c("Africa", "Australia", "Eurasia", "Papunesia", "South America", "North America"))
my_supertree = assemble_supertree(macro_groups = my_list)
