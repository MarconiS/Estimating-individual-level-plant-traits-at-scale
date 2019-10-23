test_that("using a random bag generator the pipeline creates a model with numeric AIC vector", {
  setwd("../..")
  get_random_bags(lp= 1)
  out <- pls_glm(ll = 1)
  expect_true(out$mod$AIC %>% is.numeric())
})
