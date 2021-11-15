####two predictors matrix
### parameters
n <- 100
p <- c(50,100)
nondes <- c(20,30)
sigmaondes <- c(0.05,0.02)
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data$X
X1 <- X[,(1:p[1])]
X2 <- X[,(p[1]+1):p[2]]
y <- data$y

indG <-c(rep(1,p[1]),rep(2,p[2]))

#fitting the model1
ncp <- 10
ppnu <- c(0.99,0.9)
mod.dspls <- d.spls.GLC(X=X,y=y,ncp=ncp,ppnu=ppnu,indG=indG,verbose=TRUE)
n <- dim(X)[1]
p <- dim(X)[2]

#Dimension testing
test_that("scores", { expect_equal(dim(mod.dspls$scores), c(n,ncp)) })
test_that("intercept", { expect_equal(length(mod.dspls$intercept), ncp) })
test_that("Bhat", { expect_equal(dim(mod.dspls$Bhat), c(p,ncp)) })
test_that("loadings", { expect_equal(dim(mod.dspls$loadings), c(p,ncp)) })
test_that("fitted.values", { expect_equal(dim(mod.dspls$fitted.values), c(n,ncp)) })

#residuals
test_that("residuals", { expect_setequal(mod.dspls$residuals, y-mod.dspls$fitted.values) })

#Mean of X
test_that("Xmean", { expect_setequal(apply(X, 2, mean), mod.dspls$Xmean) })

#zerovar
for (i in 2:ncp)
{
  test_that("zerovar1", { expect_gt(mod.dspls$zerovar[1,i-1],mod.dspls$zerovar[1,i]-1)})
}
test_that("zerovar1_comp1", { expect_lt(mod.dspls$zerovar[1,1], ppnu[1]*p+1)})

#zerovar
for (i in 2:ncp)
{
  test_that("zerovar2", { expect_gt(mod.dspls$zerovar[2,i-1],mod.dspls$zerovar[2,i]-1)})
}
test_that("zerovar2_comp1", { expect_lt(mod.dspls$zerovar[2,1], ppnu[2]*p+1)})

#number of variables
test_that("indG", { expect_equal(length(unique(indG)), dim(mod.dspls$zerovar)[1]) })

