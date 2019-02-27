#include "bSplineOrthoZeroMean.hpp"

// delta_x = abs(basismat[2,1] - basismat[1,1])
//     basis = basismat[,-1]
// QInt = t(basis)%*%basis*delta_x
//     bInt = colSums(basis)*delta_x
//     L = as(chol(QInt),"triangularMatrix") ## Q = L^T * L, upper tri
//     b_tilde = solve(t(L), bInt,system = "L")
//     p = dim(basis)[2]
// transMat = matrix(0,p,p-1) ##cols orthorgonal to b_tilde
//     for(j in 1:(p-1)){ 
//         transMat[j,j] = 1
//         transMat[p,j] = -b_tilde[j]/b_tilde[p]
//     }
//     transMat = qr.Q(qr(transMat))
//         transMat = solve(L, transMat,system = "L") ## transformation matrix
//         newBasisMat = basis %*% transMat
//         return(cbind(basismat[,1],newBasisMat))

orthoZeroMeanSpline::orthoZeroMeanSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_):
    bSpline( t_min_, t_max_, mOrder_, nbreaks_){

    // Compute the transformation matrix, RGlobal
    arma::mat Ddense, QInt;
    arma::vec bInt, bTilde;
    arma::mat tQ, tR, tL, tmpCoef;
    
    // Compute the integration 
    Ddense = bSpline::evalDense();
    bInt = sum(Ddense, 1) * deltaT;
    QInt = Ddense * Ddense.t() * deltaT;
    tL = arma::chol(QInt);
    
    // Create linearly independent vectors
    // orthonromal to bTilde
    bTilde = solve(tL.t(), bInt);
    tmpCoef = arma::zeros<arma::mat>(degreeOF, 
                                     degreeOF - 1);
    int k;
    for(k = 0; k < degreeOF - 1; k++){
        tmpCoef(k, k) = 1;
        tmpCoef(degreeOF - 1, k) = - bTilde(k) / bTilde(degreeOF - 1);
    }
    
    // Orthonormalized and tranform back
    arma::qr_econ(tQ, tR, tmpCoef);
    RGlobal = solve(tL, tQ);
    RGlobal = RGlobal.t(); // readly for later use
    
    Ddense = RGlobal * bSpline::evalDenseDeriv(2);
    Omega = Ddense * Ddense.t() * deltaT;
    // double nOmega = norm(Omega, "frob");
    // Omega /= nOmega;
}


// Input: evaluate spline value at t\in [t_min, t_max]
// Output: Each column is an observation
arma::mat orthoZeroMeanSpline::evalSpline(arma::vec timePoints){
    arma::mat result;
    result = RGlobal * bSpline::evalSpline(timePoints);
    return result;
}


// Output: Each column is an observation
arma::mat orthoZeroMeanSpline::evalSplineDeriv(arma::vec timePoints, int derivN){
    arma::mat result;
    result = RGlobal * bSpline::evalSplineDeriv(timePoints, derivN);
    return result;
}
