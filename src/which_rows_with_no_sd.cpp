#include <stdio.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

inline double round( double val )
{
  if( val < 0 ) return ceil(val - 0.5);
  return floor(val + 0.5);
}

double CalculateMean(double * values, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
      sum += values[i];
    }
    return (sum / size);
}

double CalculateVariance(double * values, double mean, int size)
{
    double temp = 0;

    for(int i = 0; i < size; i++)
    {
        temp += (values[i] - mean) * (values[i] - mean);
    }
    return temp / (size - 1 );
}

double Calculate_StandardDeviation(double * values, int size)
{
    double mean = CalculateMean(values, size);
    return sqrt(CalculateVariance(values, mean, size));
}

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
