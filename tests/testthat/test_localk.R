context("Local Ripley's K")

test_that("sequential local_k yields same result as spatstat::localK", {
  slk <- spatstat::localK(spatstat::cells, correction = "none", verbose = FALSE)
  lk <- local_k(spatstat::cells, correction = "none", verbose = FALSE)
  
  expect_is(slk, "fv")
  expect_is(lk, "fv")
  expect_identical(slk, lk)
})

test_that("sequential local_k with translation border correct yields same result as spatstat::localK", {
  slk <- spatstat::localK(spatstat::cells, correction = "translation", verbose = FALSE)
  lk <- local_k(spatstat::cells, correction = "translation", verbose = FALSE)
  
  expect_is(slk, "fv")
  expect_is(lk, "fv")
  expect_identical(slk, lk)
})

test_that("sequential local_k with isotropic border correct yields same result as spatstat::localK", {
  slk <- spatstat::localK(spatstat::cells, correction = "isotropic", verbose = FALSE)
  lk <- local_k(spatstat::cells, correction = "isotropic", verbose = FALSE)
  
  expect_is(slk, "fv")
  expect_is(lk, "fv")
  expect_identical(slk, lk)
})

test_that("parallel local_k with isotropic border correct yields same result as spatstat::localK", {
  cl <- parallel::makeCluster(2L)
  doParallel::registerDoParallel(cl)

  slk <- spatstat::localK(spatstat::cells, correction = "isotropic", verbose = FALSE)
  lk <- local_k(spatstat::cells, correction = "isotropic", verbose = FALSE)
  
  expect_is(slk, "fv")
  expect_is(lk, "fv")
  expect_identical(slk, lk)
  foreach::registerDoSEQ()
  cl <- parallel::stopCluster(cl)
})

test_that("local_k displays message indicating that execution is seq when parallel backend isn't registered", {
  expect_message(local_k(spatstat::cells), "No parallel backend registered. Execution will be sequential.")
})

test_that("local_k displays message indicating that progress is disabled for parallel execution", {
  cl <- parallel::makeCluster(2L)
  doParallel::registerDoParallel(cl)
  expect_message(local_k(spatstat::cells), "Progress is disabled when running in parallel.")
  foreach::registerDoSEQ()
  cl <- parallel::stopCluster(cl)
})
