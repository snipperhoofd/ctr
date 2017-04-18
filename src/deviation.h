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

#endif //__UTILITIES__
