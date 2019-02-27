#ifndef ZERO_MEAN_SPLINE_CLASS
#define ZERO_MEAN_SPLINE_CLASS


#include "bSpline.hpp"

class orthoZeroMeanSpline: public bSpline{
public:
    
    orthoZeroMeanSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_);
    
    // Input: evaluate spline value at t\in [t_min, t_max]
    // Output: Each row is an observation
    arma::mat evalSpline(arma::vec timePoints);
    arma::mat evalSplineDeriv(arma::vec timePoints, int derivN);
    
    // One less due to the zero mean constraint. 
    int getDoF(){ return degreeOF - 1; }
private:
    arma::mat RGlobal;
};

#endif

