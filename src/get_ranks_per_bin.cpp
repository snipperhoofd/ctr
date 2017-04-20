#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

std::vector<double> rankValues(std::vector<double> values, int size)
{
    std::vector<double> ranks(size);
    for(int i = 0; i < size; i++)
    {
        int rnk = 0;
        for(int z = 0; z < size; z++)
        {
            if(values[z] < values[i]){
                rnk++;
            }
        }

        ranks[i] = rnk * 1.0;
    }
    int max = *std::max_element(ranks.begin(), ranks.end());
    //divide by the max rank

    for(int i = 0; i < ranks.size(); i++)
    {
        ranks[i] = ranks[i] / max;
    }

    return ranks;
}

std::vector<double> retrieveValues(NumericVector sample,
                                std::vector<int> binContainingLines)
{
  std::vector<double> expressionValues(binContainingLines.size());
  for(int i = 0; i < binContainingLines.size(); i++)
  {
      expressionValues[i] = sample[binContainingLines[i]];
  }
  return expressionValues;
}

std::vector<int> lineNumbers(double * bins, int size, int bin)
{
  std::vector<int> lines;
  lines.reserve(size);

  for(int i = 0; i < size; i++)
  {
      if(bins[i] == bin)
      {
          lines.push_back(i);
      }
  }
  return lines;
}

NumericMatrix GetRanksPerBin(NumericMatrix samples, NumericVector bins,
                             NumericVector hqBins)
{

    NumericMatrix rankcols(samples.nrow(), samples.ncol());

    for(int i = 0; i < hqBins.size(); i++)
    {
        int bin = hqBins[i];
        std::vector<int> binContainingLines = lineNumbers(bins.begin(),
                                                          bins.size(),
                                                          bin);
        for(int j = 0; j < samples.ncol(); j++)
        {
          std::vector<double> expressionValues = retrieveValues(samples(_, j),
                                                             binContainingLines);
          std::vector<double> ranks = rankValues(expressionValues,
                                              expressionValues.size());

          for(int k = 0 ; k < binContainingLines.size(); k++)
          {
            rankcols(binContainingLines[k],j) = ranks[k];
          }

        }
    }
    return rankcols;
}



/*
 * binCOntainingLines -> sample index -> get samplevalues -> rank them
 *
 */
