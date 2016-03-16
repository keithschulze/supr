context("Utilities")

test_that("xy coords can be swapped with different coords in marks", {
  x <- c(1, 6, 8)
  y <- c(5, 8, 3)
  xc <- c(32, 76, 12)
  yc <- c(31, 35, 65)

  p <- spatstat::ppp(x, y, xrange = c(min(x), max(x)),
    yrange = c(min(y), max(y)), marks = data.frame(xc, yc))

  p_out <- exchange_xy_for_marks(p, x = "xc", y = "yc")
  expect_equal(p_out$x, xc)
  expect_equal(p_out$y, yc)
  expect_equal(p_out$marks$orig_x, x)
  expect_equal(p_out$marks$orig_y, y)
})

test_that("with_matlab execute a code block in the context 
    of an environment where the matlab is setup",
    {
        matlab_path <- "/Applications/MATLAB_R2014b.app/bin/matlab"
        expect_true(
          with_matlab(function(matlab) R.matlab::isOpen(matlab),
                      matlab_path = matlab_path)
        )
        expect_true({
          try(with_matlab(function(matlab) stop("I stopped"),
                          matlab_path = matlab_path, port=9999))
          with_matlab(function(matlab) R.matlab::isOpen(matlab),
                      matlab_path = matlab_path, port=9999)
        })
    })
