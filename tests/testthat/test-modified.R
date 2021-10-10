library(dual.spls)
library(pdist)
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data.benchmark=BCHMK(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
X <- data.benchmark$X
y <- data.benchmark$y
ncells <- 5
Datatype <- type(y,ncells=ncells)

pcal <- 70
pval <- 100-70

# nb elts in calibration for each cell
ycounts <- sapply(1:ncells,function(u) sum(Datatype==u) )
Listecal <- ceiling(ycounts*pcal/100)

index.cal <- modified.KS(X,Datatype,Listecal)

test0=which(index.cal==0)
test_that("number of calibration points", {
  expect_equal(length(index.cal), sum(Listecal))
})

test_that("content of calibration index vector", {
  expect_equal(test0, integer())
})
