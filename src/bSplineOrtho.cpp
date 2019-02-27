#include "bSplineOrtho.hpp"


orthoSpline::orthoSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_):
    bSpline( t_min_, t_max_, mOrder_, nbreaks_){

    // Compute the transformation matrix, RGlobal
    arma::mat Ddense = bSpline::evalDense();
    arma::mat QInt; 
    
    Ddense = bSpline::evalDense();
    QInt = Ddense * Ddense.t() * deltaT;
    RGlobal = arma::trans(chol(QInt)); // QInt = RGlobal^T * RGlobal

    Ddense = bSpline::evalDenseDeriv(2);
    Ddense = solve(RGlobal, Ddense);
    Omega = Ddense * Ddense.t() * deltaT;
}


// Input: evaluate spline value at t\in [t_min, t_max]
// Output: Each column is an observation
arma::mat orthoSpline::evalSpline(arma::vec timePoints){
    arma::mat result;
    result = bSpline::evalSpline(timePoints);
    result = solve(RGlobal, result);
    return result;
}


// Output: Each column is an observation
arma::mat orthoSpline::evalSplineDeriv(arma::vec timePoints, int derivN){
    arma::mat result;
    result = bSpline::evalSplineDeriv(timePoints, derivN);
    result = solve(RGlobal, result);
    return result;
    
}