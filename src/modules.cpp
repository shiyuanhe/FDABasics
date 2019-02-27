
#include "bSpline.hpp"
#include "bSplineOrtho.hpp"
#include "bSplineOrthoZeroMean.hpp"


RCPP_MODULE(MFPCA){

    
    Rcpp::class_<bSpline>("bSpline")
        .constructor<double, double, int, int>()
        .method("evalSplineDeriv", &bSpline::evalSplineDeriv)
        .method("evalSpline", &bSpline::evalSpline)
        .method("getDoF", &bSpline::getDoF)
        .method("getTMin", &bSpline::getTMin)
        .method("getTMax", &bSpline::getTMax)
        .method("getKnots", &bSpline::getKnots)
        .method("get_Omega", &bSpline::get_Omega)
    ;
    
        
    Rcpp::class_<orthoSpline>("orthoSpline")
        .derives<bSpline>("bSpline")
        .constructor<double, double, int, int>()
        .method("evalSpline", &orthoSpline::evalSpline)
        .method("evalSplineDeriv", &orthoSpline::evalSplineDeriv)
    ;
    
    Rcpp::class_<orthoZeroMeanSpline>("orthoZeroMeanSpline")
        .derives<bSpline>("bSpline")
        .constructor<double, double, int, int>()
        .method("evalSpline", &orthoZeroMeanSpline::evalSpline)
        .method("evalSplineDeriv", &orthoZeroMeanSpline::evalSplineDeriv)
        .method("getDoF", &orthoZeroMeanSpline::getDoF)
    ;
    
}
