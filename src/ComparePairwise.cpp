#include <Rcpp.h>
#include <stdio.h>
#include "statistics.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


//'@export
// [[Rcpp::export]]
List comparePairwise_C(NumericVector position_of_kegg_enzyme_A,
                       NumericVector position_of_kegg_enzyme_B,
                       NumericMatrix expressions,
                       NumericMatrix ranks)
{

    int posAlength = position_of_kegg_enzyme_A.size();
    int posBlength = position_of_kegg_enzyme_B.size();

    NumericMatrix pairwiseEuclidean(posAlength, posBlength);
    NumericMatrix pairwisePearson(posAlength, posBlength);

    for(int m = 0; m < posAlength; m++)
    {
        int pos_m = position_of_kegg_enzyme_A[m] - 1;
        NumericVector expr_A = expressions(pos_m, _);
        NumericVector rank_A = ranks(pos_m, _);
        for(int n = 0; n < posBlength; n++)
        {
            int pos_n = position_of_kegg_enzyme_B[n] - 1;
            NumericVector expr_B = expressions(pos_n, _);
            NumericVector rank_B = ranks(pos_n, _);



            double pearson = CalcNormPearson(expr_A.begin(),
                                             expr_B.begin(),
                                             expr_A.size());
            pairwisePearson(m,n) = pearson;

            pairwiseEuclidean(m,n) = CalcNormEuclidean(rank_A.begin(),
                                                       rank_B.begin(),
                                                       rank_A.size());
        }
    }
    return List::create(
      _["pairwise_correlation"] = pairwisePearson,
      _["pairwise_euclidean"] = pairwiseEuclidean
    );
}

