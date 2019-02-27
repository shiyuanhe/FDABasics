
## install GSL (install homebrew https://brew.sh/, install GSL by https://brewformulas.org/Gsl)
## install RcppGSL
## install this package

library(FDABasics)

## Spline [0,1], order = 4 (cubic spline)
## Knots 
tmin = 0
tmax = 1
order = 4
knots = 10

spline2 = new(orthoSpline, tmin, tmax, 4, 10)


tSeq = seq(tmin, tmax, length.out = 10000) ## nSample = 1000
yHat = spline2$evalSpline(tSeq)
str(yHat) ## Degree of freedom of spline X nSample

yHat[,1] #for the first sample b(t_1)\in R^{DoF}



plot(tSeq, yHat[5,], type = "l")

diag(yHat %*% t(yHat)/10000)
