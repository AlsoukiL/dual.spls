test_that("d.spls.metric works", {
  n <- 100
  p <- 50
  nondes <- 20
  sigmaondes <- 0.5
  data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)
  X <- data$X
  y <- data$y

  #fitting the model
  ncp <- 10
  ppnu<- 0.9
  mod.dspls <- d.spls.lasso(X=X,y=y,ncp=ncp,ppnu=ppnu,verbose=TRUE)
  n <- dim(X)[1]
  p <- dim(X)[2]

  predmetric= d.spls.metric(mod.dspls$fitted.values,y)

  #Dimension testing
  expect_equal(dim(mod.dspls$residuals), c(n,ncp))
  expect_equal(length(mod.dspls$fitted.values[,1]), length(y))
  expect_equal(length(predmetric$RMSE), ncp )
  expect_equal(length(predmetric$MAE), ncp)
  expect_equal(length(predmetric$Rsquared), ncp)

  i=sample(1:ncp,1)
  #RMSE
  expect_equal(
    predmetric$RMSE[i],
    sqrt(mean((y-mod.dspls$fitted.values[,i])^2)),
    tolerance = 1e-5)

  #MAE
  expect_equal(
    predmetric$MAE[i],
    mean(abs(y-mod.dspls$fitted.values[,i])),
    tolerance = 1e-5)

  #R2
  expect_equal(
    predmetric$Rsquared[i],
    1-sum((mod.dspls$residuals[,i])^2)/sum((y-mean(y))^2),
    tolerance = 1e-5)
  expect_lt(predmetric$Rsquared[i],1)

})
