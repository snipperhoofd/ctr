#include <Rcpp.h>
#include "statistics.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


class Offset {
private:
  int nrows, ncols, nmats ;

public:
  Offset( int nrows_, int ncols_, int nmats_) : nrows(nrows_),
  ncols(ncols_), nmats(nmats_){}

  int operator()( int i, int j, int k){
    return i + j * nrows + k * ( nrows * ncols ) ;
  }

} ;


int Calc_BestScoring(NumericMatrix ZscoreCorrelation,
                     NumericMatrix ZscoreEuclidean){
  double prev = 100000;//arbitrary
  int prev_idx = -1;
  int xsize = ZscoreCorrelation.ncol() * ZscoreCorrelation.nrow();
  for(int i = 0; i < xsize; i++)
  {
      double sum = ((1- ZscoreCorrelation[i]) + ZscoreEuclidean[i]);
      if(sum < prev){
        prev = sum;
        prev_idx = i;
      }
  }
  return prev_idx;
}

NumericMatrix Calc_ZscoreCorrelation(NumericMatrix correlationMatrix,
                                     List Zscores)
{
    int xsize = correlationMatrix.ncol() * correlationMatrix.nrow();
    NumericVector mu = Zscores["mu"];
    NumericVector sd = Zscores["sd"];
    for(int i = 0; i < xsize; i++)
    {
      correlationMatrix[i] = ((correlationMatrix[i] - mu[1]) / sd[1]);
    }
    return correlationMatrix;
}

NumericMatrix Calc_ZscoreEuclidean(NumericMatrix rankMatrix,
                                   List Zscores)
{
    int xsize = rankMatrix.ncol() * rankMatrix.nrow();
    NumericVector mu = Zscores["mu"];
    for(int i = 0; i < xsize; i++)
    {
        rankMatrix[i] = ((rankMatrix[i] - mu[5]) / mu[5]);
    }
    return rankMatrix;
}

std::set<int> MatchingPositions(int bin, NumericVector allBins)
{
    std::set<int> out;


    for(int i = 0; i < allBins.size(); i++)
    {
        if(allBins[i] == bin){
            out.insert(i);
        }
    }

  return out;
}
std::set<int> MatchingPositions(String KO, StringVector subsetKOS)
{
  std::set<int> out;

  for(int i = 0; i < subsetKOS.size(); i++)
  {
    if(subsetKOS[i] == KO){
      out.insert(i);
    }
  }

  return out;
}

std::vector<int> Intersection(std::set<int> a, std::set<int> b)
{
    std::set<int> intersect;
    std::set_intersection(a.begin(), a.end(),
                          b.begin(), b.end(),
                          std::inserter(intersect, intersect.begin()));
    std::vector<int> output;

    std::copy(intersect.begin(), intersect.end(), std::back_inserter(output));
    return output;
}

// [[Rcpp::export]]
List P_NRED_Distance_C(int dim_matrix,
                                StringVector subsetKOS,
                                NumericVector binNames,
                                NumericVector allBins,
                                StringVector allKOs,
                                NumericMatrix expression,
                                NumericMatrix ranks,
                                List Z_scores) {




    Offset offset(dim_matrix, dim_matrix, subsetKOS.size());

    NumericVector PairwiseBinArrayPearson = NumericVector(Dimension(dim_matrix,
                                                                    dim_matrix,
                                                                    subsetKOS.size()));

    NumericVector PairwiseBinArrayEuclidean = NumericVector(Dimension(dim_matrix,
                                                               dim_matrix,
                                                               subsetKOS.size()));

    //Positions in the genome where the scores come from
    NumericVector Pairwise_PositionsA = NumericVector(Dimension(dim_matrix,
                                                               dim_matrix,
                                                               subsetKOS.size()));
    NumericVector Pairwise_PositionsB = NumericVector(Dimension(dim_matrix,
                                                               dim_matrix,
                                                               subsetKOS.size()));

    NumericVector v1_exp(6);
    NumericVector v2_exp(6);
    for(int x = 0; x < dim_matrix; x++)
    {
      printf("x");
        std::set<int> positionsGenomeA = MatchingPositions(binNames[x],
                                                         allBins);
        for(int y = 0; y < dim_matrix; y++)
        {
            std::set<int> positionsGenomeB = MatchingPositions(binNames[y],
                                                             allBins);
            for(int z = 0; z <subsetKOS.size(); z++)
            {
              std::set<int> allPositionsKegg = MatchingPositions(
                                                          subsetKOS[z], allKOs);
              std::vector<int> positionKeggA = Intersection(allPositionsKegg,
                                                    positionsGenomeA);
              std::vector<int> positionKeggB = Intersection(allPositionsKegg,
                                                    positionsGenomeB);

              if (positionKeggA.size() != 0 && positionKeggB.size() != 0){
                NumericMatrix maxPairwiseCorrelation(positionKeggA.size(),
                                                      positionKeggB.size());
                NumericMatrix maxPairwiseEuclidean(positionKeggA.size(),
                                                      positionKeggB.size());

                for(int m = 0; m < positionKeggA.size(); m++)
                {
                  int positionKeggA_current = positionKeggA[m];
                  for(int n = 0; n < positionKeggB.size(); n++)
                  {
                    int positionKeggB_current = positionKeggB[n];
                    v1_exp = expression(positionKeggA_current, _);
                    v2_exp = expression(positionKeggB_current, _);
                    maxPairwiseCorrelation(m,n) = CalcNormPearson(v1_exp.begin(),
                                                 v2_exp.begin(),
                                                 v1_exp.size());

                    NumericVector v1_rank = ranks(positionKeggA_current, _);
                    NumericVector v2_rank = ranks(positionKeggB_current, _);
                    maxPairwiseEuclidean(m,n) = CalcNormEuclidean(v1_rank.begin(),
                                                   v2_rank.begin(),
                                                   v1_rank.size());
                  }
                }


                //convert Z scores
                NumericMatrix ZscoresCorrelation = Calc_ZscoreCorrelation(maxPairwiseCorrelation,
                                                                          Z_scores);
                NumericMatrix ZscoreEuclidean = Calc_ZscoreEuclidean(maxPairwiseEuclidean,
                                                                     Z_scores);
                int bestScoring = Calc_BestScoring(ZscoresCorrelation, ZscoreEuclidean);

                PairwiseBinArrayPearson[offset(x,y,z)] = maxPairwiseCorrelation[bestScoring];
                PairwiseBinArrayEuclidean[offset(x,y,z)] = maxPairwiseEuclidean[bestScoring];



              }
          }
        }
    }

   // NumericMatrix Zscore_pairwise_gene_correlation = Calc_ZscoreCorrelation(PairwiseBinArrayPearson,
    //                                                                        Z_scores);


  return List::create(
   _["Pairwise Euclidean"] = PairwiseBinArrayEuclidean,
   _["Pairwise Pearson"] = PairwiseBinArrayPearson);
//   _['pearsons'] = PairwiseBinArrayPearson,
//   _['nred'] = ZscoreEuclidean
// );
}


