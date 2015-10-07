context("Spatial analysis")

test_that("all_types should use specified marks as type factors", {
  points <- read.moleculelist_nstorm("./data/nstormfull.txt")
  win008a <- spatstat::owin(xrange = c(27000, 29000), yrange = c(18000, 20000),
    unitname = "nm")
  points_crop <- points[win008a]
  lc <- all_types(points_crop, fun="Lcross", which.marks = "channel.name", 
    correction = "isotropic", envelope = TRUE)
  expect_is(lc, "fasp")
  expect_error(all_types(points_crop, fun="Lcross",
    which.marks = c("channel.name", "xc"), correction = "isotropic",
    envelope = TRUE), "which.marks should be a string specifying a single column of marks.")
})


test_that("multi_local_l with mark.original=TRUE returns a ppp marked with 
  L(r) value at each coord", {
  cl <- parallel::makeCluster(2L)
  doParallel::registerDoParallel(cl)

  points_llc <- multi_local_l(spatstat::amacrine, correction = "isotropic",
    rvalue = 0.15, mark.original=TRUE, verbose=FALSE)
  expect_is(points_llc, "ppp")
  expect_true("L0.15" %in% colnames(spatstat::marks(points_llc)))
  
  points_llc <- multi_local_l(points_llc, which.marks = "spatstat::amacrine",
    correction = "isotropic", rvalue = 0.05, mark.original=TRUE, verbose=FALSE)
  expect_true("L0.05" %in% colnames(spatstat::marks(points_llc)))

  foreach::registerDoSEQ()
  cl <- parallel::stopCluster(cl)
})

test_that("multi_local_lcross with mark.original=TRUE returns a ppp marked with 
  L(r)cross values at each coord", {
  cl <- parallel::makeCluster(2L)
  doParallel::registerDoParallel(cl)

  points_llc <- multi_local_lcross(spatstat::amacrine, correction = "isotropic",
    rvalue = 0.15, mark.original=TRUE, verbose=FALSE)
  expect_is(points_llc, "ppp")
  expect_true("L0.15" %in% colnames(spatstat::marks(points_llc)))
  foreach::registerDoSEQ()
  cl <- parallel::stopCluster(cl)
})

test_that("co_cluster_l_cross with mark.original=TRUE returns a ppp marked with 
  L(r) and L(r) cross values at each coord", {
  cl <- parallel::makeCluster(2L)
  doParallel::registerDoParallel(cl)

  points_llc <- co_cluster_l_cross(spatstat::amacrine, correction = "isotropic",
    rvalue = 0.15, mark.original=TRUE, verbose=FALSE)
  expect_is(points_llc, "ppp")
  expect_true("L0.15" %in% colnames(spatstat::marks(points_llc)))
  expect_true("L0.15c" %in% colnames(spatstat::marks(points_llc)))

  points_llc <- co_cluster_l_cross(points_llc, which.marks = "spatstat::amacrine",
    correction = "isotropic", rvalue = 0.05, mark.original=TRUE, verbose=FALSE)
  expect_true("L0.05" %in% colnames(spatstat::marks(points_llc)))
  expect_true("L0.05c" %in% colnames(spatstat::marks(points_llc)))

  foreach::registerDoSEQ()
  cl <- parallel::stopCluster(cl)
})
