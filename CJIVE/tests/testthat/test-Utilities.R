test_that("MatVar works", {
  X = matrix(rnorm(9),3,3)
  expect_equal(MatVar(X), sum(diag(t(X)%*%X)))
})

test_that("MatVar2 works", {
  X = matrix(rnorm(9),3,3)
  expect_equal(MatVar2(X), sum(X^2))
})


test_that("Chordal norm works", {
  X = matrix(rnorm(15),5,3)
  X2 = X%*%rbind(diag(rnorm(2)),rep(0,2))
  expect_equal(chord.norm.diff(X,X2), 0)
})
