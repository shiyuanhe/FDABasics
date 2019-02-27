#include "bSpline.hpp"

bSpline::bSpline(double t_min_, double t_max_, int mOrder_, int nbreaks_){
    mOrder = mOrder_;
    t_min = t_min_;
    t_max = t_max_;
    nbreaks = nbreaks_;
    degreeOF = mOrder + nbreaks - 2;
    
    // For dense evaluation and integration approximation.
    seqN = 20000;
    deltaT =  (t_max - t_min) / static_cast<double>(seqN);
    
    if(!bw) gsl_bspline_free(bw);
    bw = gsl_bspline_alloc(mOrder, nbreaks);
    gsl_bspline_knots_uniform(t_min, t_max, bw);

    // Get knots
    knots = arma::vec(bw->knots->size, arma::fill::zeros);

    for(int i = 0; i < knots.size(); i++){
        knots[i] = gsl_vector_get(bw->knots, i);
    }
    
}



// Input: evaluate spline value at t\in [t_min, t_max]
// Output: Each column is an observation
arma::mat bSpline::evalSpline(arma::vec timePoints){
    arma::mat result(degreeOF, timePoints.n_elem);
    for(int i=0; i < timePoints.n_elem; i++){
        result.col(i) = bSpline::evalSplineSingle(timePoints(i));
    }
    return result;
}


// Output: Each column is an observation
arma::mat bSpline::evalSplineDeriv(arma::vec timePoints, int derivN){
    arma::mat result(degreeOF, timePoints.n_elem);
    for(int i=0; i < timePoints.n_elem; i++){
        result.col(i) = bSpline::evalSplineDerivSingle(timePoints(i), derivN);
    }
    return result;
    
}


// evaluate spline value at a single t\in [t_min, t_max]
arma::vec bSpline::evalSplineSingle(double t){
    RcppGSL::vector<double> result_tmp(degreeOF);
    arma::vec result(degreeOF);
    gsl_bspline_eval(t, result_tmp, bw);
    Rcpp::NumericVector tmp(result_tmp.begin(), result_tmp.end());
    result = arma::vec(tmp);
    return result;
}


// evaluate spline derivative value at t\in [t_min, t_max]
arma::vec bSpline::evalSplineDerivSingle(double t, int derivN){
    RcppGSL::matrix<double> result_tmp(degreeOF, derivN + 1);
    gsl_bspline_deriv_eval(t, derivN, result_tmp, bw);
    arma::vec result(degreeOF);
    for(int i = 0; i < degreeOF; i++){
        result(i) = result_tmp(i, derivN);//lastCol.vector.data[i];
    }
    return result;
}


// evaluate over a dense sequence of T
arma::mat bSpline::evalDense(){
    arma::mat result(degreeOF, seqN);
    arma::vec seqT = arma::linspace<arma::vec>(t_min, t_max, seqN);
    for(int i =0; i<seqN; i++){
        result.col(i) = evalSplineSingle(seqT(i));
    }
    return result;
}


// evaluate over a dense sequence of T
arma::mat bSpline::evalDenseDeriv(int derivN){
    arma::mat result(degreeOF, seqN);
    arma::vec seqT = arma::linspace<arma::vec>(t_min, t_max, seqN);
    for(int i =0; i<seqN; i++){
        result.col(i) = evalSplineDerivSingle(seqT(i), derivN);
    }
    return result;
}
