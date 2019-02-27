# FDABasics


This package is the basics of functional data analysis.

1. B-spline, orthonormal B-spline, orthonormal B-spline with zero-mean?
2. Fit the mean function (or surface).


## Example

The example construct 4th order B-spline on the interval [tmin, tmax], with 10 equally spaced knots. 

```r
library(FDABasics)

# parameters
tmin = 0
tmax = 1
order = 4
nknots = 10

## construct spline
spline1 = new(bSpline, tmin,  tmax, order, nknots)

## spline 
nSeq = 100
tSeq = seq(tmin, tmax, length.out = nSeq)
y1 = spline1$evalSpline(tSeq)
plot(tSeq, y1[1,], type = "l")

# get the degree of freedom
spline1$getDoF()
```

The evaluated spline values are stored into the columns of `y1`.

> <span style="color:red">
The spline evaluation output is column wise.
 Each column corresponds to an observation. 
</span>


Construct orthonormal B-spline.

```r
spline2 = new(orthoSpline, tmin, tmax, order, nknots)
```
Construct orthonormal B-spline with zero mean.

```r
spline3 = new(orthoZeroMeanSpline, tmin, tmax, oroder, nknots)
```

