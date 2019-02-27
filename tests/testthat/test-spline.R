context("Test Three Splines")
spline1 = new(bSpline, 0, 1, 4, 10)
spline2 = new(orthoSpline, 0, 1, 4, 10)
spline3 = new(orthoZeroMeanSpline, 0, 1, 4, 10)
nSeq = 20000
tSeq = seq(0, 1, length.out = nSeq)
y1 = spline1$evalSpline(tSeq)

test_that("correct dimension",{
    K = spline2$getDoF()
    Km1 = spline3$getDoF()
    expect_equal(K - 1, Km1)
})

test_that("Orthonormal Spline",{
    y2 = spline2$evalSpline(tSeq)
    Z = diag(rep(1,spline2$getDoF()))
    QInt2 = y2 %*% t(y2) / nSeq
    err = sum(abs(Z - QInt2))
    expect_lt(err, 1e-5)
})

test_that("Orthonormal Spline with zero mean",{
    y3 = spline3$evalSpline(tSeq)
    Z = diag(rep(1, spline3$getDoF()))
    QInt3 = y3 %*% t(y3) / nSeq
    err = sum(abs(Z - QInt3))
    expect_lt(err, 1e-5)
    
    bInt3 = rowSums(y3) / nSeq
    expect_lt(sum(abs(bInt3)), 1e-5)
})



