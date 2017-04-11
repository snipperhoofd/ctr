#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
#include <stdlib.h>
#include <vector>

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


class StdDeviation
{

private:
    int max;
    double value[100];
    double mean;



public:
    double CalculateMean()
    {
        double sum = 0;
        for(int i = 0; i < max; i++)
        {
            sum += value[i];
        }
        return (sum / max);
    }



    double CalculateVariane()
    {
        mean = CalculateMean();
        double temp = 0;
        for(int i = 0; i < max; i++)
        {
             temp += (value[i] - mean) * (value[i] - mean) ;
        }
        return temp / max;
    }


    double CalculateSampleVariane()
    {
        mean = CalculateMean();
        double temp = 0;
        for(int i = 0; i < max; i++)
        {
             temp += (value[i] - mean) * (value[i] - mean) ;
        }
        return temp / (max - 1);
    }



    int SetValues(double *p, int count)
    {
        if(count > 100)
            return -1;
        max = count;
        for(int i = 0; i < count; i++)
            value[i] = p[i];
        return 0;
    }


    double Calculate_StandardDeviation()
    {
        return sqrt(CalculateVariane());
    }



    double Calculate_SampleStandardDeviation()
    {
        return sqrt(CalculateSampleVariane());
    }

};

class Calculator
{

    private:
    double XSeries[100];
    double YSeries[100];
    int max;

    StdDeviation x;
    StdDeviation y;



    public:
    void SetValues(double *xvalues, double *yvalues, int count)
    {
        for(int i = 0; i < count; i++)
        {
            XSeries[i] = xvalues[i];
            YSeries[i] = yvalues[i];
        }
        x.SetValues(xvalues, count);
        y.SetValues(yvalues, count);
        max = count;
    }




    double Calculate_Covariance()
    {
        double xmean = x.CalculateMean();
        double ymean = y.CalculateMean();
        double total = 0;

        for(int i = 0; i < max; i++)
        {
            total += (XSeries[i] - xmean) * (YSeries[i] - ymean);
        }
        return total / max;
    }


    double Calculate_Correlation()
    {
        double cov = Calculate_Covariance();
        double correlation = cov / (x.Calculate_StandardDeviation() * y.Calculate_StandardDeviation());
        return correlation;
    }

};


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
