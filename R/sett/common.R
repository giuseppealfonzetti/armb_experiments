sims_setup <- tidyr::expand_grid(
  n = c(5e3, 1e4, 5e4, 1e5),
  p = c(50, 100, 200, 500),
  seed = 1:5
  )

qs::qsave(sims_setup, file = "output/sims_setup.qs")

