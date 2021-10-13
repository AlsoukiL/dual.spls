library(pdist)

#Data1
vec=1:100

#Positive value
test_that("norm1", {expect_gt(norm1(vec),0)})
test_that("norm2", {expect_gt(norm2(vec),0)})
#equal norm1 and norm2

test_that("norm1", {expect_equal(norm1(vec),sum(abs(vec)))})
test_that("norm2", {expect_equal(norm2(vec),sqrt(sum(vec^2)))})




#Data2
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
X <- data.benchmark$X
y <- data.benchmark$y

#Positive value
test_that("norm1", {expect_gt(norm1(y),0)})
test_that("norm2", {expect_gt(norm2(y),0)})
#equal norm1 and norm2

test_that("norm1", {expect_equal(norm1(y),sum(abs(y)))})
test_that("norm2", {expect_equal(norm2(y),sqrt(sum(y^2)))})


