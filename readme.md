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

## Orthonormal B-spline

Evaluate with Cholesky decompostion
$$
\mathbf{Q} = \int \mathbf{b}(t) \mathbf{b}(t) ^Tdt
 = \mathbf{L}^T\mathbf{L}.
$$
Suppose $\mathbf{e}_k$ is the standard Euclidean basis, then
$$
\boldsymbol{\beta}_k = \mathbf{L}^{-1}\mathbf{e}_k.
$$
The spline coef is
$$
\mathbf{B} = (\beta_1, \cdots, \beta_{K})
= \mathbf{L}^{-1}\cdot \mathbf{I},
$$
and
$$
\mathbf{B}^T\mathbf{b}(t) = \mathbf{I}\cdot(\mathbf{L}^T)^{-1} \mathbf{b}(t)
$$

<span style="color:red">
**Issue** The evaluated $\Omega$ has very large entries, and it may have negative eigenvalues.
</span>

## Orthonormal B-spline with zero-mean

We want to find the spline coef
$$
\mathbf{B} = (\beta_1, \cdots, \beta_{K-1})
$$
such that
$$
\begin{align}
\mathbf{0} =&\int \mathbf{b}^T(t)\boldsymbol{\beta}_i dt   = \mathbf{b}^T\boldsymbol{\beta}_i\\
\delta_{ij} =& \int 
\mathbf{b}^T(t)\boldsymbol{\beta}_i
\mathbf{b}^T(t)\boldsymbol{\beta}_j dt
= \boldsymbol{\beta}_i^T\mathbf{Q}
\boldsymbol{\beta}_j
\end{align}
$$
where
$\int \mathbf{b}(t) dt = \mathbf{b}$
and 
$\mathbf{Q} = \int \mathbf{b}(t) \mathbf{b}(t) ^Tdt$.

Suppose we have the Cholesky decompostion
$\mathbf{Q} = \mathbf{L}^T\mathbf{L}$, 
and take the tranformation
$$
\begin{align}
\tilde{\mathbf{b}} & = (\mathbf{L}^T)^{-1} \mathbf{b}\\
\tilde{\boldsymbol{\beta}}_i 
&= \mathbf{L}\boldsymbol{\beta}_i.
\end{align}
$$
Then we are trying to solve the system
$$
\tilde{\mathbf{b}}^T\tilde{\boldsymbol{\beta}}_i  = 0 
\text{ and }
\tilde{\boldsymbol{\beta}}_i^T
\tilde{\boldsymbol{\beta}}_j = \delta_{ij}.
$$
In the code, we find a group of linearly independent vectors orthonormal to 
$\tilde{\mathbf{b}}$, 
$$
\left(0,\cdots,\ \underset{\text{at k}}{1}\ ,0,\cdots, ´
-\frac{\tilde{b}_k}{\tilde{b}_K}
\right)^T
$$
then make them orthonormal to each other by QR decomposition. Solve the equation to get the coef in the original scale
$$
\mathbf{L}^{-1}\tilde{\boldsymbol{\beta}}_i 
= \boldsymbol{\beta}_i
$$


