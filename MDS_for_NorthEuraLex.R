# set working directory
setwd("C:/Users/stein/Dropbox/Studium/7. Semester/BA-Thesis/BA-Thesis_NorthEuraLex")

# load libraries
# library(gstat)
# library(sf)
# library(spData)
# library(terra)
library(tidyverse)
# library(tmap)

conceptdata_loc = "northeuralex-0.9-forms.tsv"
geodata_loc = "northeuralex-0.9-language-data.tsv"

if(!file.exists(conceptdata_loc)) {
  download.file(
    "http://www.sfs.uni-tuebingen.de/~jdellert/northeuralex/0.9/northeuralex-0.9-forms.tsv",
    dest = conceptdata_loc
  )
}

if(!file.exists(geodata_loc)) {
  download.file(
    "http://www.sfs.uni-tuebingen.de/~jdellert/northeuralex/0.9/northeuralex-0.9-language-data.tsv",
    dest = geodata_loc
  )
}

# read in the NorthEuraLex data
conceptdata = read_tsv("northeuralex-0.9-forms.tsv")
geodata = read_tsv("northeuralex-0.9-language-data.tsv")

conceptdata = 
  conceptdata %>% 
  # select the 3 relevant columns
  select(Language_ID, Concept_ID, ASJP) %>%
  # add the missing row
  add_row(Language_ID = "xal", Concept_ID = "Auge::N", ASJP = "nyudEN")

# make a vector with the ASJP word list
asjp_concepts = c("Auge::N", "Ohr::N", "Nase::N", "Zahn::N", "Zunge::N",
                  "Busen::N", "Hand::N", "Knie::N", "Haut::N", "Blut::N", 
                  "Knochen::N", "Leber::N", "Sonne::N", "Stern::N", "Wasser::N",
                  "Stein::N", "Feuer::N", "Berg::N", "Baum::N", "Blatt::N",
                  "Horn::N", "Hund::N", "Fisch::N", "Laus::N", "Mensch::N", 
                  "Name::N", "Pfad::N", "Nacht::N", "voll::A", "neu::A", 
                  "ich::PRN", "du::PRN", "wir::PRN", "eins::NUM", "zwei::NUM",
                  "trinken::V", "sterben::V", "kommen::V", "sehen::V", 
                  "hören::V")

# filter out other concepts
conceptdata = filter(conceptdata, Concept_ID %in% asjp_concepts)

# get all the languages 
languages = 
  geodata %>% select(iso_code)
#%>% filter(family == "Indo-European") 

# for each language, compile a list of the concepts and the corresponding word(s) for that concept
language_params = list()
for (i in 1:nrow(languages)) {
  lang = 
    conceptdata %>%
    filter(Language_ID == languages$iso_code[i]) %>%
    select(Concept_ID, ASJP)
  language_params[[i]] = lang
}


# create a matrix to store the pairwise levenshtein distances
levenshtein_distances = 
  matrix(1, nrow(languages), nrow(languages)) # add-1 smoothing


# iterate through all pairs of languages
for (i in 1:(nrow(languages) - 1)) {
  # print to see the program is actually doing something 
  # (because the loop takes a while)
  print(i)
  for (j in (i + 1):nrow(languages)) {
    # creates a tibble comparing two languages
    join = inner_join(language_params[[i]], language_params[[j]], 
                      by=c("Concept_ID" = "Concept_ID"))
    language_distance = 0
    k = 1
    n_concepts = 40
    remove_concept = TRUE
    for (concept in asjp_concepts) {
      n_words = 0
      concept_distance = 0
      while ((join$Concept_ID[k] == concept) && k != nrow(join) + 1) {
        remove_concept = FALSE
        n_words = n_words + 1
        concept_distance = 
          concept_distance + 
          (adist(join$ASJP.x[k], join$ASJP.y[k]) /
             max(nchar(join$ASJP.x[k]), nchar(join$ASJP.y[k])))
        k = k + 1
      }
      if (remove_concept) {n_concepts = n_concepts - 1} 
      else {
        remove_concept = TRUE
        language_distance = language_distance + (concept_distance / n_words)
      }
      counter = 0
      concept_distance = 0
    }
    levenshtein_distances[i, j] = (language_distance / n_concepts) 
    levenshtein_distances[j, i] = (language_distance / n_concepts)
  }
}


write.csv(levenshtein_distances, "levenshtein_distances.csv", row.names = FALSE)


pmiDists = levenshtein_distances

hc <- pmiDists %>%
  as.dist() %>%
  hclust()

ggdendrogram(hc, rotate = T, size=.05)


nClasses = 9

soundClasses = cutree(hc, k=nClasses) %>% as.numeric()
l = languages %>% .[[1]]


for (i in 1:nClasses) {
  print(paste0(l[soundClasses == i], collapse=" "))
}




# apply multi-dimensional scaling on the resulting matrix
multi_dimensional_scaling = cmdscale(levenshtein_distances, k=3) %>%
  as_tibble()

mds = 
  as_tibble(multi_dimensional_scaling) %>%
  mutate(languages = languages$iso_code) %>%
  mutate(x = V1 - min(.$V1), y = V2 - min(.$V2), z = V3 - min(.$V3)) %>%
  mutate(r = x / max(.$x), g = y / max(.$y), b = z / max(.$z)) %>%
  mutate(col = rgb(r,g,b)) %>%
  select(languages, col, r, g, b)


geodata = geodata %>% 
  rename(languages = iso_code) %>%
  select(languages, latitude, longitude) %>%
  inner_join(mds)

world %>%
  ggplot() +
  geom_sf() +
  geom_point(data=geodata, aes(x=longitude, y = latitude), col=geodata$col, size=1.5)

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

# remove everything on the 60? meridian from our world data, effectively cutting all polygons that cross the meridian
world_without_55 = 
  world %>% 
  st_difference(meridian_55_west) %>% 
  st_transform("+proj=eqearth lon_0=125")
  

# convert the coordinates in geodata into simple features
geodata_sf = 
  geodata %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>%
  st_transform("+proj=eqearth lon_0=125")

geodata_sf %>%
  tm_shape() +
  tm_symbols(size=0.2, col="col", border.lwd=1, alpha=1) +
  tm_shape(world_without_55) +
  tm_polygons(alpha=0, border.alpha = 0.3)