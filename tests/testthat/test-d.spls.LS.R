test_that("d.spls.LS works", {

  data=d.spls.NIR
  X=data$NIR
  X=as.matrix(X)
  y=data$density
  n <- length(y)
  p <- dim(X)[2]

  # fitting the model
  ncp <- 5
  ppnu <- 0.9
  expect_warning(d.spls.LS(X=X,y=y,ncp=ncp,ppnu=ppnu,verbose=TRUE),"deflated XtX is close to being singular on component number 5")
  ncp=4
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
  expect_equal(apply(X, 2, mean), mod.dspls$Xmean)

  # zerovar
  for (i in 2:ncp)
  {
    expect_gt(mod.dspls$zerovar[i-1],mod.dspls$zerovar[i]-1)
  }
  expect_lt(mod.dspls$zerovar[1], ppnu*p+1)})

