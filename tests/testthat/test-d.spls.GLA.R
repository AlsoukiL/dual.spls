test_that("d.spls.GLA works", {

  #### two predictors matrix
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

  # fitting the model
  ncp <- 10
  ppnu <- c(0.99,0.9)
  mod.dspls <- d.spls.GLA(X=X,y=y,ncp=ncp,ppnu=ppnu,indG=indG,verbose=TRUE)
  n <- dim(X)[1]
  p <- dim(X)[2]

  # dimension testing
  expect_equal(dim(mod.dspls$scores), c(n,ncp))
  expect_equal(length(mod.dspls$intercept), ncp)
  expect_equal(dim(mod.dspls$Bhat), c(p,ncp))
  expect_equal(dim(mod.dspls$loadings), c(p,ncp))
  expect_equal(dim(mod.dspls$fitted.values), c(n,ncp))

  # residuals
  expect_equal(mod.dspls$residuals, y-mod.dspls$fitted.values, tolerance = 1e-5)

  # mean of X
  expect_equal(apply(X, 2, mean), mod.dspls$Xmean, tolerance = 1e-5)

  # zerovar
  for (i in 2:ncp)
  {
    expect_gt(mod.dspls$zerovar[1,i-1],mod.dspls$zerovar[1,i]-1)
  }
  expect_lt(mod.dspls$zerovar[1,1], ppnu[1]*p+1)

  # number of variables
  expect_equal(length(unique(indG)), dim(mod.dspls$zerovar)[1])

})

