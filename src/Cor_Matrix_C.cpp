#include <Rcpp.h>
#include <math.h>
#include "statistics.h"
using namespace Rcpp;


std::vector<double> subtract(double * rank_v1, double * rank_v2, int size)
{
    std::vector<double> subtracted_list(size);
    for(int i = 0; i < size; i++)
    {
        subtracted_list[i] = (rank_v1[i] - rank_v2[i]);
    }
    return subtracted_list;
}

double sqrt_sum_multiply(std::vector<double> s, int size)
{
    double sum = 0.0;
    for(int i = 0; i < size; i++)
    {
        sum += s[i] * s[i];
    }
    return sqrt(sum);
}

//'@export
// [[Rcpp::export]]
List Cor_Matrix_C(NumericVector v1, NumericVector v2,
                           NumericMatrix expression, NumericMatrix ranks) {
    NumericMatrix pairwise_gene_correlation(v1.size(), v2.size());
    NumericMatrix pairwise_gene_euclidean(v1.size(), v2.size());
    for(int i = 0; i < v1.size(); i++)
    {
        for(int j = 0; j < v2.size(); j++)
        {
            NumericVector sub_v1 = expression(v1[i], _);
            NumericVector sub_v2 = expression(v2[j], _);
            pairwise_gene_correlation(i, j) = CalcNormPearson(sub_v1.begin(),
                                                              sub_v2.begin(),
                                                              sub_v1.size());

            NumericVector rank_v1 = ranks(v1[i],_ );
            NumericVector rank_v2 = ranks(v2[j],_ );

            std::vector<double> subtract_list = subtract(rank_v1.begin(), rank_v2.begin(), rank_v1.size());
            pairwise_gene_euclidean(i, j ) =  sqrt_sum_multiply(subtract_list, subtract_list.size());
        }
    }

    return List::create(
      _["correlation"] = pairwise_gene_correlation,
      _["euclidean"] = pairwise_gene_euclidean
    );
}

