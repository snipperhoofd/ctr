#ifndef __UTILITIES__
#define __UTILITIES__

inline double CalculateMean(double * values, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
      sum += (values[i] * 1.0);
    }
    return (sum / size);
}


inline double CalculateVariance(double * values, double mean, int size)
{
    double temp = 0;

    for(int i = 0; i < size; i++)
    {
        temp += ((values[i] * 1.0) - mean) * ((values[i] * 1.0) - mean);
    }
    return temp / (size - 1 );
}


inline double Calculate_StandardDeviation(double * values, int size)
{
    double mean = CalculateMean(values, size);
    return sqrt(CalculateVariance(values, mean, size));
}

inline double CalcNormPearson(double * genomeA, double * genomeB, int n)
{
  double r = 0.0;

  double xbar = CalculateMean(genomeA, n);
  double ybar = CalculateMean(genomeB, n);
  double xstdev = Calculate_StandardDeviation(genomeA, n);
  double ystdev = Calculate_StandardDeviation(genomeB, n);

  for(int i = 0; i < n; i++)
  {
    r += (((genomeA[i] - xbar) / xstdev) * ((genomeB[i] - ybar) / ystdev));
  }
  return (r /= (n-1));
}

inline double CalcNormEuclidean(double * vector1, double * vector2, int size)
{
  double total = 0.0;
  for(int i = 0; i < size; i++)
  {
    total +=  ((vector1[i]-vector2[i]) * (vector1[i] - vector2[i]));
  }
  return sqrt(total);
}


#endif //__UTILITIES__
