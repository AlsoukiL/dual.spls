
n <- 200
p <- 100
nondes <- 150
sigmaondes <- 0.01
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data$X
y <- data$y

#fitting the model
ncp=10
ppnu=0.9
mod.dspls <- d.spls.ridge(X=X,y=y,ncp=ncp,ppnu=ppnu,nu2=0.05,verbose=TRUE)
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
  test_that("zerovar", { expect_gt(mod.dspls$zerovar[i-1],mod.dspls$zerovar[i]-1)})
}
test_that("zerovar1", { expect_lt(mod.dspls$zerovar[1], ppnu*p+1)})

