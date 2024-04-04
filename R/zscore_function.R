
# Calculate z-scores for all simulated variables, set direction of treatment
# effect, generate average z-score for each observation over all sims

zscore_function <- function(datalist,
                            direction_vector) {

  datalist.combined <- lapply(datalist, function(df) {
    data.new <- df %>%
      dplyr::select(starts_with("V")) %>%
      mutate(across(everything(), ~(scale(.) %>% as.vector))) %>%
      rename_with(~paste0("z.", .), starts_with("V")) %>%
      sweep(2, FUN = "*", direction_vector) %>%
      mutate(z.average = rowMeans(.))

    cbind.data.frame(df, data.new)
  })
  return(datalist.combined)
}
