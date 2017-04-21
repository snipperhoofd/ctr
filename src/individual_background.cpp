#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include "statistics.h"

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
 * Selects the row counts of genes in a bin
 *
 *
 */
std::vector<int> positionsOfGenome(StringVector allBins, char bin);


bool compareKO(String a, String b)
{
    if( a != "" && a == b)
      return true;
    else
      return false;
}

List getMatchingKO(std::vector<int> positionsOfGenomeA,
                                       std::vector<int> positionsOfGenomeB,
                                       StringVector KOTerms)
{

  std::unordered_map<std::string, std::unordered_set<int> > mapA;
  std::unordered_map<std::string, std::unordered_set<int> > mapB;
  std::unordered_set<int> test;


  for(int i = 0; i < positionsOfGenomeA.size(); i++){
    String koA = KOTerms[positionsOfGenomeA[i]];
    for(int j = 0; j < positionsOfGenomeB.size(); j++)
    {
      if(compareKO(koA, KOTerms[positionsOfGenomeB[j]]))
      {
        if(mapA.find(koA) != mapA.end())
        {
          mapA[koA].insert(positionsOfGenomeA[i]);
          mapB[koA].insert(positionsOfGenomeB[j]);
        }
        else{
          std::unordered_set<int> setA;
          setA.reserve(positionsOfGenomeA.size());
          std::unordered_set<int> setB;
          setB.reserve(positionsOfGenomeB.size());

          setA.insert(positionsOfGenomeA[i]);
          setB.insert(positionsOfGenomeB[j]);
          mapA[koA] = setA;
          mapB[koA] = setB;
        }
      }
    }
  }



  return List::create(
    _["mapA"] = mapA,
    _["mapB"] = mapB
  );
}






List sample_KO(std::vector<int> positionsOfGenomeA,
                                   std::vector<int> positionsOfGenomeB,
               StringVector KOTerms)
{
  List values = getMatchingKO(positionsOfGenomeA,
                                                 positionsOfGenomeB,
                                                 KOTerms);


 // std::string randomKO = matchingTerms[rand() % matchingTerms.size()];
  //which lines do have that random term in A and B

  //take a random shared KO
  return List::create(
    _["sample_KO_position_of_A"] = NULL,
    _["sample_KO_position_of_B"] = NULL,
    _["KO_positions_of_A"] = NULL,
    _["KO_positions_of_B"] = NULL
  );
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
        //Take two random genomes
        int genomeABinIndex = rand() % HighQBins.size();
        int genomeBBinIndex = rand() % HighQBins.size();

        //determine the positions of the genes associated with a certain genome/bin
        std::vector<int> positionsOfGenomeA = positionsOfGenome(allBins.begin(),
                                                                allBins.size(),
                                                                HighQBins[genomeABinIndex]);
        std::vector<int> positionsOfGenomeB = positionsOfGenome(allBins.begin(),
                                                                allBins.size(),
                                                                HighQBins[genomeBBinIndex]);
        //within the randome genomes take a random gene each
        int randomGeneA = rand() % positionsOfGenomeA.size();
        int randomGeneB = rand() % positionsOfGenomeB.size();

        //Pearson
        NumericVector geneA = RNAseqExpressionCounts(randomGeneA, _);
        NumericVector geneB = RNAseqExpressionCounts(randomGeneB, _);
        RandomPairwiseGenePearson[i] = CalcNormPearson(geneA.begin(), geneB.begin(),
                                                       geneA.size());

        //NRED
        geneA = RNAseqExpressionRanks(randomGeneA, _);
        geneB = RNAseqExpressionRanks(randomGeneB, _);
        RandomPairwiseGeneEuclidean[i] = CalcNormEuclidean(geneA.begin(),
                                                           geneB.begin(),
                                                           geneA.size());

        List test = sample_KO(positionsOfGenomeA, positionsOfGenomeB,  KOTerms);

    }

    return List::create(
     _["random_pairwise_gene_pearson"] = RandomPairwiseGenePearson,
     _["random_pairwise_gene_euclidean"] = RandomPairwiseGeneEuclidean
   );
}
