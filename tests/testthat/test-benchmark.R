library(dual.spls)
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
X <- data.benchmark$X
y <- data.benchmark$y
y0 <- data.benchmark$y0

test_that("right dimensions for X", {
  expect_equal(dim(X)[1], n)
  expect_equal(dim(X)[2], p)
})

test_that("right dimensions for y", {
  expect_equal(length(y), n)
  expect_equal(length(y), length(y0))
})

test_that("right parameters", {
  expect_equal(data.benchmark$sigmaondes, sigmaondes)
})
