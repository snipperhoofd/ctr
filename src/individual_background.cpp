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


std::vector<int> positionsOfGenome(double * allBins, int size, int bin)
{
    std::vector<int> positions;
    positions.reserve(20000);//arbitrary value


    for(int i = 0; i < size; i++){
        if(allBins[i] == bin){
            positions.push_back(i);
        }
    }
    return positions;
}

double CalcNormPearson(NumericVector genomeA, NumericVector genomeB)
{
    int n = genomeA.size(); //should be equal to genomeB size
    double r = 0.0;

    double xbar = CalculateMean(genomeA.begin(), n);
    double ybar = CalculateMean(genomeB.begin(), n);
    double xstdev = Calculate_StandardDeviation(genomeA.begin(), n);
    double ystdev = Calculate_StandardDeviation(genomeB.begin(), n);

    for(int i = 0; i < n; i++)
    {
      r += (((genomeA[i] - xbar) / xstdev) * ((genomeB[i] - ybar) / ystdev));
    }
    return (r /= (n-1));
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
                              NumericVector allBins,
                              NumericVector HighQBins,
                              int N)
{
    NumericVector RandomPairwiseGenePearson(N);
    NumericVector RandomPairwiseGeneEuclidean(N);

    for(int i = 0; i < N; i++)
    {
        int genomeABinIndex = rand() % HighQBins.size();
        int genomeBBinIndex = rand() % HighQBins.size();

        //determine the positions of the genes associated with a certain genome/bin
        std::vector<int> positionsOfGenomeA = positionsOfGenome(allBins.begin(),
                                                                allBins.size(),
                                                                HighQBins[genomeABinIndex]);
        std::vector<int> positionsOfGenomeB = positionsOfGenome(allBins.begin(),
                                                                allBins.size(),
                                                                HighQBins[genomeBBinIndex]);

        int randomGeneA = rand() % positionsOfGenomeA.size();
        int randomGeneB = rand() % positionsOfGenomeB.size();

        //correlation
        NumericVector geneA = RNAseqExpressionCounts(randomGeneA, _);
        NumericVector geneB = RNAseqExpressionCounts(randomGeneB, _);

        RandomPairwiseGenePearson[i] = CalcNormPearson(geneA, geneB);
    }

    return List::create(
     _["random_pairwise_gene_pearson"] = RandomPairwiseGenePearson,
     _["random_pairwise_gene_euclidean"] = RandomPairwiseGeneEuclidean
   );
}
