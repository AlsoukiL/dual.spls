test_that("d.spls.predict works", {

### parameters
n <- 100
p <- 50
nondes <- 20
sigmaondes <- 0.5
data=d.spls.simulate(n=n,p=p,nondes=nondes,sigmaondes=sigmaondes)

X <- data$X
y <- data$y

pcal <- 70
ncells <- 3

split <- d.spls.calval(X=X,pcal=pcal,y=y,ncells=ncells)

indcal= split$indcal
indval= split$indval

Xcal=X[indcal,]
Xval=X[indval,]
ycal=y[indcal]
yval=y[indval]

#fitting the model
ncp=10
mod.dspls <- d.spls.lasso(X=Xcal,y=ycal,ncp=ncp,ppnu=0.9,verbose=TRUE)

ycalhat=mod.dspls$fitted.values
rescal=mod.dspls$residuals
# predictions on validation
liste_cp=1:ncp

yvalhat=d.spls.predict(mod.dspls,Xval, liste_cp=liste_cp)

#Dimension testing
expect_equal(dim(yvalhat), c(length(indval),length(liste_cp)))

#Prediction testing
expect_equal(
  ycalhat, d.spls.predict(mod.dspls,Xcal, liste_cp=1:ncp),
  tolerance = 1e-5
)

})

