test_that("d.spls.calval works", {

n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
X <- data$X
y <- data$y
ncells <- 5
Datatype <- d.spls.type(y,ncells=ncells)

pcal <- 70

# nb elts in calibration for each cell
ycounts <- sapply(1:ncells,function(u) sum(Datatype==u) )
Listecal <- ceiling(ycounts*pcal/100)

index.cal <- d.spls.calval(X,Datatype,Listecal)

test0=which(index.cal==0)

# number of calibration points
expect_equal(length(index.cal), sum(Listecal))


# content of calibration index vector"
  expect_equal(test0, integer())
})
