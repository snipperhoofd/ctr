#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
#include <stdlib.h>
#include "deviation.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


/**
 * Calculates the Pearson Correlation between two vectors
 *
 * @param vector1 vector containing doubles
 * @param vector2 vector containing doubles
 * @param size Integer representing the size of both vectors (equal size)
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
double CalcNormPearson(double * vector1, double * vector2, int size);


/**
 * Calculates the Normalized Euclidean distance of two vectors.
 *
 * @param vector1 vector containing doubles
 * @param vector2 vector containing doubles
 * @param size Integer representing the size of both vectors (equal size)
 * @return the normalised euclidean distance
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
double CalcNormEuclidean(double * vector1, double * vector2, int size);

/**
 * Determines if two samples share the same KO term
 *
 * @param K01 char array containing the first KO term
 * @param K02 char array containing the second KO term
 * @return true if the terms match
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
bool SampleKO(StringVector * KOterms);

/**
 * Selects the row counts of genes in a bin
 *
 *
 */
std::vector<int> positionsOfGenome(StringVector allBins, char bin);


bool SampleKO(String a, String b)
{
    if(a != "" && a == b)
      return true;
    else
      return false;
}


std::vector<int> positionsOfGenome(StringVector allBins, char bin)
{

    std::vector<int> positions;
    positions.reserve(20000);//arbitrary value


    for(int i = 0; i < allBins.size(); i++){
        if(allBins[i] == bin){
            positions.push_back(i);
        }
    }
    return positions;
}

double CalcNormPearson(NumericVector genomeA, NumericVector genomeB)
{
  /*
  Calculator calc;


    int size = genomeA.size();


    double x[size];
    double y[size];

    //Transform NumericVector to an array of doubles
    //std::vector<double> x = Rcpp::as<std::vector<double>>(genomeA);
    //std::vector<double> y = Rcpp::as<std::vector<double>>(genomeB);
    //Make the calculator accept std::vector<double> would be good
    for(int i = 0; i < size;i++){
        x[i] = genomeA[i];
        y[i] = genomeB[i];
    }

*/
  //calc.SetValues(x, y, size);
  //double correlation = calc.Calculate_Correlation();
  return 1.0;
}

double CalcNormEuclidean(double * vector1, double * vector2, int size)
{
  int total = 0;
    for(int i = 0; i < size; i++)
    {
        total +=  ((vector1[i]-vector2[i]) * (vector1[i] - vector2[i]));
    }
    return sqrt(total);
}

/*todo:
 * KO_pairwise_gene_pearson
 * KO_pairwise_gene_euclidean
 * H_KO_pairwise_gene_pearson
 * H_KO_pairwise_gene_euclidean
 * H_random_pairwise_gene_pearson
 * H_random_pairwise_gene_euclidean
 */

//' @title Individual_KO_background
//' @export
// [[Rcpp::export]]
List Individual_KO_background(NumericMatrix RNAseqExpressionCounts,
                              NumericMatrix RNAseqExpressionRanks,
                              StringVector KOTerms,
                              StringVector allBins,
                              StringVector HighQBins,
                              int N)
{

    NumericVector RandomPairwiseGenePearson(N);
    NumericVector RandomPairwiseGeneEuclidean(N);

    for(int i = 0; i < N; i++)
    {
        int genomeABinIndex = rand() % HighQBins.size();
        int genomeBBinIndex = rand() % HighQBins.size();

        //determine the positions of the genes associated with a certain genome/bin
        std::vector<int> positionsOfGenomeA = positionsOfGenome(allBins, HighQBins[genomeABinIndex][0]);
        std::vector<int> positionsOfGenomeB = positionsOfGenome(allBins, HighQBins[genomeBBinIndex][0]);

        int randomGeneA = rand() % positionsOfGenomeA.size();
        int randomGeneB = rand() % positionsOfGenomeB.size();


        printf("value: %d", positionsOfGenomeA[randomGeneA]);
      //  NumericVector gA = RNAseqExpressionCounts(positionsOfGenomeA[randomGeneA], _);
      //  NumericVector gB = RNAseqExpressionCounts(positionsOfGenomeB[randomGeneB], _);
      //  RandomPairwiseGenePearson[i] = CalcNormPearson(RNAseqExpressionCounts(positionsOfGenomeA[randomGeneA], _),
        //                                               RNAseqExpressionCounts(positionsOfGenomeB[randomGeneB], _));
       // RandomPairwiseGeneEuclidean[i] = CalcNormEuclidean(x, y, genomeA.size());

    }



    return List::create(
     _["random_pairwise_gene_pearson"] = RandomPairwiseGenePearson,
     _["random_pairwise_gene_euclidean"] = RandomPairwiseGeneEuclidean
   );
}
