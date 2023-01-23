library(gstat)
library(sf)
library(spData)
library(terra)
library(tidyverse)
library(tmap)
library(glottoTrees)

setwd("C:/Users/stein/Dropbox/Studium/7. Semester/BA-Thesis/BA-Thesis_NorthEuraLex")

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
conceptdata_raw = read_tsv("northeuralex-0.9-forms.tsv")
geodata = read_tsv("northeuralex-0.9-language-data.tsv")

# exclude rows which are still under review
conceptdata = conceptdata_raw %>% filter(Next_Step == "validate")

conceptdata = 
  conceptdata %>% 
  # select the 3 relevant columns
  select(Language_ID, Concept_ID, ASJP)

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

# for each language, compile a list of the concepts and the corresponding word(s) for that concept
language_params = list()
for (i in 1:nrow(languages)) {
  lang = 
    conceptdata %>%
    filter(Language_ID == languages$iso_code[i]) %>%
    select(Concept_ID, ASJP)
  language_params[[i]] = lang
}


PMI_data_loc = "pnas.1500331112.sd04.csv"

if(!file.exists(PMI_data_loc)) {
  download.file(
    "http://www.pnas.org/lookup/suppl/doi:10.1073/pnas.1500331112/-/DCSupplemental/pnas.1500331112.sd04.csv",
    dest = PMI_data_loc
  )
}


PMI_scores = read.table("pnas.1500331112.sd04.csv", sep = ",", check.names=FALSE)
head(PMI_scores)

asjp_to_index = c("!" = 1, "3" = 2, "4" = 3, "5" = 4, "7" = 5, "8" = 6, "C" = 7,
                  "E" = 8, "G" = 9, "L" = 10, "N" = 11, "S" = 12, "T" = 13, 
                  "X" = 14, "Z" = 15, "a" = 16, "b" = 17, "c" = 18, "d" = 19, 
                  "e" = 20, "f" = 21, "g" = 22, "h" =23, "i" = 24, "j" = 25, 
                  "k" = 26, "l" = 27, "m" = 28, "n" = 29, "o" = 30, "p" = 31, 
                  "q" = 32, "r" = 33, "s" = 34, "t" = 35, "u" = 36, "v" = 37, 
                  "w" = 38, "x" = 39, "y" = 40, "z" = 41)

get_PMI_score <- function(c1, c2){
  PMI_scores[[asjp_to_index[c1], asjp_to_index[c2]]]
}

affine = function(word1, word2, similarity_func = get_PMI_score, 
                  gap_start = -2.4930, gap_continuation = -1.7057) {
  
  length1 = nchar(word1)
  length2 = nchar(word2)
  
  m = matrix(0, length1 + 1, length2 + 1)
  x = matrix(0, length1 + 1, length2 + 1)
  y = matrix(0, length1 + 1, length2 + 1)
  
  for (i in 2:(length1 + 1)){
    m[i, 1] = -Inf
    x[i, 1] = gap_start + (i - 1) * gap_continuation
    y[i, 1] = -Inf
  } 
  
  for (j in 2:(length2 + 1)){
    m[1, j] = -Inf
    x[1, j] = -Inf
    y[1, j] = gap_start + (j - 1) * gap_continuation
  }
  
  for (i in 2:(length1 + 1)){
    for (j in 2:(length2 + 1)){
      m[i,j] = similarity_func(substr(word1, i - 1, i - 1), substr(word2, j - 1, j - 1)) +
        max(m[i - 1, j - 1], x[i - 1, j - 1], y[i - 1, j - 1])
      
      x[i, j] = max(gap_start + m[i - 1, j], gap_continuation + x[i - 1, j])
      
      y[i, j] = max(gap_start + m[i, j - 1], gap_continuation + y[i, j - 1])      
    }
  }
  
  max(m[length1 + 1, length2 + 1], x[length1 + 1, length2 + 1], y[length1 + 1, length2 + 1])
}

affine_gap_distances = 
  matrix(0, nrow(languages), nrow(languages))


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
          ((affine(join$ASJP.x[k], join$ASJP.y[k])) 
           / max(nchar(join$ASJP.x[k]), nchar(join$ASJP.y[k])))
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
    affine_gap_distances[i, j] = (language_distance / n_concepts) 
    affine_gap_distances[j, i] = (language_distance / n_concepts)
  }
}

write.csv(affine_gap_distances, "language_distances_averaged.csv", row.names = FALSE)

