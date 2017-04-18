/**
 * Calculate the mean of an array with doubles
 *
 * @param values A pointer to an array of doubles
 * @param size The length of the array
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
double CalculateMean(double * values, int size)
{
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
      sum += values[i];
    }
    return (sum / size);
}

/**
 * Generic function to calculate the sample variance
 *
 * @param values A pointer to an array of doubles
 * @param mean The mean value of the array with doubles
 * @param size The lenght of the array
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
double CalculateVariance(double * values, double mean, int size)
{
    double temp = 0;

    for(int i = 0; i < size; i++)
    {
        temp += (values[i] - mean) * (values[i] - mean);
    }
    return temp / (size - 1 );
}

/**
 * Generic function to calculate the standard deviation as the square root of
 * the sample variance.
 *
 * @param values A pointer to an array of doubles
 * @param size The length of the array
 *
 * @author Joris van Steenbrugge
 * @version 1.0 10 april 2017
 */
double Calculate_StandardDeviation(double * values, int size)
{
    double mean = CalculateMean(values, size);
    return sqrt(CalculateVariance(values, mean, size));
}
