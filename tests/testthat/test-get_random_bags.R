test_that("testing get_random_bags ... expecting one pixel per tree crown", {
  setwd("../..")
  random_bags <- get_random_bags()
  dims <- dim(random_bags)[1]
  expect_gte(dims, 1)
})
