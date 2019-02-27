#ifndef BSPLINE_CLASS
#define BSPLINE_CLASS


// Generate ordinary   B-Spline

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <cmath>

class bSpline{
public:
    virtual ~bSpline(){
        if(!bw) gsl_bspline_free(bw);
    }
    
    bSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_);
    
    
    // evaluate over a dense sequence of T
    // Output: Each column is an observation
    virtual arma::mat evalSpline(arma::vec timePoints);
    virtual arma::mat evalSplineDeriv(arma::vec timePoints, int derivN);
    
    
    // Get parameters for R
    virtual int getDoF(){ return degreeOF; }
    double getTMin(){ return t_min; }
    double getTMax(){ return t_max; }
    arma::vec getKnots(){ return knots; }
    arma::mat get_Omega(){ return Omega; }
    
    
protected:
    //seqN, deltaT for dense evaluation and integration approximation.
    int mOrder, degreeOF, nbreaks, seqN;
    double t_min, t_max, deltaT;
    arma::mat Omega;
    arma::vec knots;
    gsl_bspline_workspace *bw;

    // on a dense sequence of t on [t_min, t_max]
    arma::mat evalDense();
    arma::mat evalDenseDeriv(int derivN);
    // evaluate spline derivative value at t\in [t_min, t_max]
    arma::vec evalSplineDerivSingle(double t, int derivN);
    // evaluate spline value at t\in [t_min, t_max]
    arma::vec evalSplineSingle(double t);
    
};


#endif

