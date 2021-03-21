

double slopeLimiter_MC(double *arrayTemp, int j)
{
    //declare variables
    double minmod_a;
    double MC_ratio_L;
    double minmod_b;
    double minmod_b_L;
    double minmod_a_L;
    double sigma;
    double sigma_L;
    double MC_ratio;
    double min_a;
    double min_b;
    double min_a_L;
    double min_b_L;
    
    //Monotonized Central Slope Limiter
    MC_ratio = (arrayTemp[j] - arrayTemp[j-1])/(arrayTemp[j+1] - arrayTemp[j]);

    //if (2*MC_ratio < 2 && 2*MC_ratio < 0.5*(1+MC_ratio))

    min_a = fmin(2*MC_ratio, 0.5*(1+MC_ratio));
    min_b = fmin(0.5*(1+MC_ratio),2);

    sigma = fmax(0, fmin(min_a,min_b));

    MC_ratio_L = (arrayTemp[j-1] - arrayTemp[j-2])/(arrayTemp[j] - arrayTemp[j-1]);

    min_a_L = fmin(2*MC_ratio_L, 2);
    min_b_L = fmin(2*MC_ratio_L, 0.5*(1+MC_ratio_L));

    sigma_L = fmax(0, fmin(min_a_L, min_b_L));

    return sigma - sigma_L;
}