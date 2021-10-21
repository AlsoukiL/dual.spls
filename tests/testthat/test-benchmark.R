library(dual.spls)

####one predictors matrix
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

####two predictors matrix
### parameters
n <- 100
p <- c(50,100)
nondes <- c(20,30)
sigmaondes <- c(0.05,0.02)
data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data.benchmark$X
X1 <- X[,(1:p[1])]
X2 <- X[,(p[1]+1):p[2]]
y <- data.benchmark$y

test_that("right dimensions for X", {
  expect_equal(dim(X)[1], n)
  expect_equal(dim(X)[2], sum(p))
})

test_that("right dimensions for y", {
  expect_equal(length(y), n)
  expect_equal(length(y), length(y0))
})

test_that("right parameters", {
  expect_equal(data.benchmark$sigmaondes, sigmaondes)
})

test_that("warning test",{
  expect_warning(BCHMK(n=n,p=c(100,200),nondes=50,sigmaondes=c(0.005,0.5)))
  })
