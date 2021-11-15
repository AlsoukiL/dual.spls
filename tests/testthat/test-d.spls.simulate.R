library(dual.spls)

####one predictors matrix
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data1=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
Xa <- data1$X
ya <- data1$y
y0a <- data1$y0

test_that("right dimensions for Xa", {
  expect_equal(dim(Xa)[1], n)
  expect_equal(dim(Xa)[2], p)
})

test_that("right dimensions for ya", {
  expect_equal(length(ya), n)
  expect_equal(length(ya), length(y0a))
})

test_that("right parameters", {
  expect_equal(data1$sigmaondes, sigmaondes)
})

####two predictors matrix
### parameters
n <- 100
p <- c(50,100)
nondes <- c(20,30)
sigmaondes <- c(0.05,0.02)
data2=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

Xb <- data2$X
X1 <- Xb[,(1:p[1])]
X2 <- Xb[,(p[1]+1):p[2]]
yb <- data2$y
y0b <- data2$y0

test_that("right dimensions for Xb", {
  expect_equal(dim(Xb)[1], n)
  expect_equal(dim(Xb)[2], sum(p))
})

test_that("right dimensions for yb", {
  expect_equal(length(yb), n)
  expect_equal(length(yb), length(y0b))
})

test_that("right parameters", {
  expect_equal(data2$sigmaondes, sigmaondes)
})

test_that("warning test",{
  expect_warning(d.spls.simulate(n=n,p=c(100,200),nondes=50,sigmaondes=c(0.005,0.5)))
  })
