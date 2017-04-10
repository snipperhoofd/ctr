#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
#include <stdlib.h>

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
 * @return A char array containg
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
char* SharedKO(char * KO1, char * KO2);


class StdDeviation{

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

class Calculator{

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




double CalcNormPearson(double * vector1, double * vector2, int size)
{
  Calculator calc;
  calc.SetValues(vector1, vector2, size);
  return calc.Calculate_Correlation();
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
 * I was trying to make a shared KO function but it might be over complicating
 * KO_pairwise_gene_pearson
 * KO_pairwise_gene_euclidean
 * H_KO_pairwise_gene_pearson
 * H_KO_pairwise_gene_euclidean
 * H_random_pairwise_gene_pearson
 * H_random_pairwise_gene_euclidean
 */

// [[Rcpp::export]]
List Individual_KO_background(NumericMatrix RNAseqExpressionCounts,
                              NumericMatrix RNAseqExpressionRanks,
                              StringVector KOTerms,
                              int N)
{

    NumericVector RandomPairwiseGenePearson(N);
    NumericVector RandomPairwiseGeneEuclidean(N);

    for(int i = 0; i < N; i++)
    {
        int indexA = rand() % N;
        int indexB = rand() % N;

        NumericVector genomeA = RNAseqExpressionCounts(indexA, _);
        NumericVector genomeB = RNAseqExpressionCounts(indexB, _);

        double x[genomeA.size()];
        double y[genomeB.size()];

        //Transform NumericVector to an array of doubles
        for(int i = 0; i < genomeA.size();i++){
            x[i] = genomeA[i];
            y[i] = genomeB[i];
        }



        RandomPairwiseGenePearson[i] = CalcNormPearson(x,y, genomeA.size());
        RandomPairwiseGeneEuclidean[i] = CalcNormEuclidean(x, y, genomeA.size());

        //Check if random KO terms are equal
        if(KOTerms[indexA] == KOTerms[indexA]){
          printf("TRUE")
        }
    }

   return List::create(
     _["random_pairwise_gene_pearson"] = RandomPairwiseGenePearson,
     _["random_pairwise_gene_euclidean"] = RandomPairwiseGeneEuclidean
   );
}
