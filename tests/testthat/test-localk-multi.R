context("Local multitype Ripley's K")

test_that("the mean value of cross-type local K values should be the same as Kcross for a given r", {
  kc <- spatstat::Kcross(spatstat::amacrine, correction="isotropic", r=c(0, 0.05, 0.1, 0.15))
  lkc0.05 <- local_k_cross(spatstat::amacrine, correction="isotropic", rvalue=0.05, verbose=FALSE)
  lkc0.1 <- local_k_cross(spatstat::amacrine, correction="isotropic", rvalue=0.1, verbose=FALSE)
  lkc0.15 <- local_k_cross(spatstat::amacrine, correction="isotropic", rvalue=0.15, verbose=FALSE)

  expect_equal(mean(lkc0.05), kc$iso[2])
  expect_equal(mean(lkc0.1), kc$iso[3])
  expect_equal(mean(lkc0.15), kc$iso[4])
})
