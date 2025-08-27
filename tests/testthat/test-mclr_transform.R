test_that("mclr.transform", {
  Z <- matrix(c(1, 2, 0, 4, 5, 6, 0, 8, 9), nrow = 3, byrow = TRUE)
  transformed_Z <- mclr.transform(Z, base = 10, eps = 0.1)
  expect_true(nrow(transformed_Z) == 3)
  expect_true(ncol(transformed_Z) == 3)
  expect_true(is.matrix(transformed_Z))
})
