context("Visualisation")

# ggplot suport for spatstat
test_that("single type ppp unpacks in ggplot", {
  plt <- ggplot2::ggplot(spatstat::cells)

  expect_that(plt, is_a("ggplot"))
  expect_equal(ncol(plt$data), 2)
  expect_equal(plt$data, as.data.frame(spatstat::cells))
})

test_that("multi-type ppp unpacks in ggplot", {
  plt <- ggplot2::ggplot(spatstat::lansing)

  expect_is(plt, "ggplot")
  expect_equal(ncol(plt$data), 3)
  expect_equal(colnames(plt$data), c("x", "y", "marks"))
  expect_true(length(levels(plt$data$marks)) == 6)
})

test_that("multi-variable multi-type ppp unpacks in ggplot", {
  plt <- ggplot2::ggplot(spatstat::gorillas)

  expect_equal(ncol(plt$data), 5)
  expect_equal(colnames(plt$data), c("x", "y", "group", "season", "date"))
})

test_that("psp objects unpack in ggplot", {
  plt <- ggplot2::ggplot(spatstat::copper$SouthLines)

  expect_is(plt, "ggplot")
  expect_equal(ncol(plt$data), 4)
  expect_equal(colnames(plt$data), c("x0", "y0", "x1", "y1"))
})

test_that("fv objects unpack in ggplot", {
  k_lansing <- spatstat::Kest(spatstat::lansing)
  expect_is(k_lansing, "fv")

  plt <- ggplot2::ggplot(k_lansing)
  expect_is(plt, "ggplot")
  expect_equal(colnames(plt$data)[1], "r")
  expect_equal(colnames(plt$data)[2], "theo")
})

test_that("im objects unpack in ggplot", {
  im <- density(spatstat::lansing)

  plt <- ggplot2::ggplot(im)
  expect_is(plt, "ggplot")
  expect_equal(colnames(plt$data), c("x", "y", "value"))
})

test_that("rectangular owin objects unpack in ggplot", {
  win <- spatstat::owin(xrange = c(0, 1000), yrange = c(0, 1000))

  plt <- ggplot2::ggplot(win)
  expect_equal(colnames(plt$data), c("xmin", "xmax", "ymin", "ymax"))
  expect_equal(nrow(plt$data), 1)
})
