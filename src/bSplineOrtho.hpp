#ifndef ORTHO_SPLINE_CLASS
#define ORTHO_SPLINE_CLASS


#include "bSpline.hpp"

class orthoSpline: public bSpline{
public:
    orthoSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_);
    
    // Input: evaluate spline value at t\in [t_min, t_max]
    // Output: Each row is an observation
    arma::mat evalSpline(arma::vec timePoints);
    arma::mat evalSplineDeriv(arma::vec timePoints, int derivN);
    
private:
    arma::mat RGlobal;
};

#endif

