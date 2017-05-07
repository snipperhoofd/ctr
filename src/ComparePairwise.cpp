#include <Rcpp.h>
#include "statistics.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


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
        NumericVector expr_A = expressions(position_of_kegg_enzyme_A[m], _);
        NumericVector rank_A = ranks(position_of_kegg_enzyme_A[m], _);
        for(int n = 0; n < posBlength; n++)
        {
            NumericVector expr_B = expressions(position_of_kegg_enzyme_A[n], _);
            NumericVector rank_B = ranks(position_of_kegg_enzyme_A[n], _);

            pairwisePearson(m,n) = CalcNormPearson(expr_A.begin(),
                                                   expr_B.begin(),
                                                   expr_A.size());
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

