
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
#Error testing
test_that("Error testing", {expect_error(calval(X,pcal=pcal,ncells=10))})
test_that("Error testing", {expect_error(calval(X,Datatype=Datatype,ncells=10))})

L1=floor(0.7*length(which(Datatype==1)))
L2=floor(0.8*length(which(Datatype==2)))
L3=floor(0.6*length(which(Datatype==3)))
L4=floor(0.7*length(which(Datatype==4)))
L5=floor(0.8*length(which(Datatype==5)))

Listecal=c(L1,L2,L3,L4)

expect_error(calval(X,Datatype=Datatype,Listecal=Listecal))

calval=calval(X,y=y,pcal=pcal,Datatype=Datatype,ncells=5)
index.cal=calval$indcal
test0=which(index.cal==0)

# nb elts in calibration for each cell
ycounts=sapply(1:ncells,function(u) sum(Datatype==u) )
Listecal=ceiling(ycounts*pcal/100)

test_that("number of calibration points", {
  expect_equal(length(index.cal), sum(Listecal))
})

test_that("content of calibration index vector", {
  expect_equal(test0, integer())
})

