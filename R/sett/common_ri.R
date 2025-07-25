sims_setup <- tidyr::expand_grid(
  n = c(100, 200, 500),
  p = c(5, 25, 50),
  sz = c(20),
  seed = 1
) |>
  dplyr::mutate(
    m = 10
  )

qs::qsave(sims_setup, file = "output/sims_ri_setup.qs")
