test_that("d.spls.LS works", {

  set.seed(14)
  n <- 1000
  p <- 40
  nondes <- 100
  sigmaondes <- 0.01
  data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
  X <- data$X
  y <- data$y

  # fitting the model
  ncp <- 4
  ppnu <- 0.9
  expect_warning(d.spls.LS(X=X,y=y,ncp=ncp,ppnu=ppnu,verbose=TRUE),"deflated XtX is close to being singular on component number 4")
  ncp=3
  mod.dspls <- d.spls.LS(X=X,y=y,ncp=ncp,ppnu=ppnu,verbose=TRUE)
  n <- dim(X)[1]
  p <- dim(X)[2]

  # dimension testing
  expect_equal(dim(mod.dspls$scores), c(n,ncp))
  expect_equal(length(mod.dspls$intercept), ncp)
  expect_equal(dim(mod.dspls$Bhat), c(p,ncp))
  expect_equal(dim(mod.dspls$loadings), c(p,ncp))
  expect_equal(dim(mod.dspls$fitted.values), c(n,ncp))

  # residuals
  expect_setequal(mod.dspls$residuals, y-mod.dspls$fitted.values)

  # mean of X
  expect_setequal(apply(X, 2, mean), mod.dspls$Xmean)

  # zerovar
  for (i in 2:ncp)
  {
    expect_gt(mod.dspls$zerovar[i-1],mod.dspls$zerovar[i]-1)
  }
  expect_lt(mod.dspls$zerovar[1], ppnu*p+1)})

