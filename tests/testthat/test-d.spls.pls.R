test_that("d.spls.pls works", {

  n <- 100
  p <- 50
  nondes <- 20
  sigmaondes <- 0.5
  data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
  X <- data$X
  y <- data$y

  # fitting the model
  ncp <- 10
  mod.dspls <- d.spls.pls(X=X,y=y,ncp=ncp,verbose=TRUE)
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

})



