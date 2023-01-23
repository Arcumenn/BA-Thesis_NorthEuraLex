pairwise_matrices = list()

affine_v = Vectorize(affine)

# create a matrix to store the pairwise levenshtein distances
avg_levenshtein_dists = 
  matrix(1, nrow(languages), nrow(languages)) # add-1 smoothing

for (i in 1:(2 - 1)) {
  # print to see the program is actually doing something 
  # (because the loop takes a while)
  print(i)
  for (j in (i + 1):2) {
    # creates a tibble comparing two languages
    joined_tibble = 
      inner_join(language_params[[i]], language_params[[j]], 
                      by=c("Concept_ID" = "Concept_ID")) %>%
      mutate(avg_levenshtein_dist = adist_v(ASJP.x, ASJP.y) / 
                                    pmax(nchar(ASJP.x), nchar(ASJP.y))) %>%
      group_by(Concept_ID) %>%
      slice_min(order_by = avg_levenshtein_dist, with_ties = F)
      
    # avg_levenshtein_dists[i, j] = avg_levenshtein_dists[j, i] = 
        # mean(joined_tibble$avg_levenshtein_dist)
    
    gamma = 0
    for (k in 1:(nrow(joined_tibble) - 1)) {
      for (l in (k + 1):nrow(joined_tibble)) {
        gamma = 
          gamma + 
          adist(joined_tibble$ASJP.x[k], joined_tibble$ASJP.y[l]) / 
          max(nchar(joined_tibble$ASJP.x[k]), nchar(joined_tibble$ASJP.y[l]))
      }
    }
    gamma = gamma / (nrow(joined_tibble) * (nrow(joined_tibble) - 1) / 2)
    avg_levenshtein_dists[i, j] = avg_levenshtein_dists[j, i] =
      mean(joined_tibble$avg_levenshtein_dist) / gamma
  }
}

write.csv(avg_levenshtein_dists, "ldn.csv", row.names = FALSE)


dissimilarities = matrix(1, nrow(languages), nrow(languages))

for (i in 1:(nrow(languages) - 1)) {
  print(i)
  for (j in (i + 1):nrow(languages)) {
    joined_tibble = 
      inner_join(language_params[[i]], language_params[[j]], by=c("Concept_ID" = "Concept_ID")) %>%
      mutate(similarity = affine_v(ASJP.x, ASJP.y) / 
               pmax(nchar(ASJP.x), nchar(ASJP.y))) %>%
      group_by(Concept_ID) %>%
      slice_min(order_by = similarity, with_ties = F)
    
    m = matrix(0, nrow(joined_tibble), nrow(joined_tibble))
    diag(m) = NaN
    for (k in 1:nrow(joined_tibble)) {
      for (l in 1:nrow(joined_tibble)) {
        if(k == l){next}
        m[k, l] = affine(joined_tibble$ASJP.x[k], joined_tibble$ASJP.y[l]) / 
          max(nchar(joined_tibble$ASJP.x[k]), nchar(joined_tibble$ASJP.y[l]))
      }
    }
    
    non_synonymous_wordpairs = nrow(joined_tibble) * (nrow(joined_tibble) - 1)
    for (k in 1:nrow(joined_tibble)){
      numerator = 1 + sum(m > joined_tibble$similarity[k], na.rm = T)
      # print(numerator)
      denominator = 1 + non_synonymous_wordpairs
      # print(denominator)
      calibrated_sim = -log(numerator / denominator)
      # print(calibrated_sim)
      joined_tibble$similarity[k] = calibrated_sim
    }
    dissimilarities[i, j] = dissimilarities[j, i] = 
      log(log(non_synonymous_wordpairs)) - log(mean(joined_tibble$similarity))
  }
}
