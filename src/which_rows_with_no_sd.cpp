#include <stdio.h>
#include <Rcpp.h>
#include <math.h>
#include "deviation.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

/**
 * Round a floating point number by adding or substracting 0.5 before rounding.
 *
 * @param val A double containing the value
 *
 * @author Joris van Steenbrugge
 * @version 1.0 13 april 2017
 */
inline double round( double val )
{
  if( val < 0 ) return ceil(val - 0.5);
  return floor(val + 0.5);
}



//'@export
// [[Rcpp::export]]
NumericVector which_rows_with_no_sd_cpp(NumericMatrix x)
{
    NumericVector out(x.nrow());

    for(int i = 0; i < x.nrow(); i++)
    {
        NumericVector expressions(x.ncol());
        expressions = x(i, _);
        double stdev = Calculate_StandardDeviation(expressions.begin(), x.ncol());
        stdev = (round(stdev * 1000.0) / 1000.0);
        out(i) = stdev;
    }
    return out;
}
