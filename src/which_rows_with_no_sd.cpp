#include <stdio.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]



double CalculateMean(double * values, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
      sum += values[i];
    }
    return (sum / size);
}

double CalculateVariance(double * values, double mean, int size){
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
StringMatrix which_rows_with_no_sd_cpp(StringMatrix x, NumericVector sampleCols)
{
    int colSize = sampleCols.size();
    StringMatrix out(x.nrow(), x.ncol());




    for(int i = 0; i < x.nrow(); i++)
    {
        double expressions[colSize];
        //Extract an array containing the sample expression values
        for(int j = 0; j < colSize; j++)
        {
            std::string a = as<std::string>(x(i,sampleCols[j]));
            double expression = atof(a.c_str());
            expressions[j] = expression;
        }

        int stdev = Calculate_StandardDeviation(expressions, colSize) + 0.5;
        if( stdev != 0)
        {
          out(i, _) = x(i, _);
        }

    }
    return out;
}
