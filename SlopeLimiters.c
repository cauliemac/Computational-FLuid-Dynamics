#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include"functions.h"

//Monotonized Central slope limiter
double slopeLimiter_MC(double *arrayTemp, int j)
{
    //declare variables
    double sigma;
    double sigma_L;
    double MC_ratio;
    double MC_ratio_L;
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

//Minmod slope limiter
double slopeLimiter_minmod(double *arrayTemp, int j)
{   
    /*
    a is j-1
    b is j
    c is j+1
    d is j-2
    */

    double minmod_a;
    double minmod_b;
    double minmod_b_L;
    double minmod_a_L;
    double sigma;
    double sigma_L;
    double dx = 2/500;

    minmod_a = (arrayTemp[j] - arrayTemp[j-1])/dx;
    minmod_b = (arrayTemp[j+1] - arrayTemp[j])/dx;

    if(minmod_a * minmod_b < 0)
        sigma = 0;
    else if (abs(minmod_a) < abs(minmod_b) && minmod_a * minmod_b > 0)
        sigma = minmod_a;
    else if (abs(minmod_b) < abs(minmod_a) && minmod_a * minmod_b > 0)
        sigma = minmod_b;

    minmod_a_L = (arrayTemp[j-1] - arrayTemp[j-2])/dx;
    minmod_b_L = (arrayTemp[j] - arrayTemp[j-1])/dx;

    if(minmod_a_L * minmod_b_L < 0)
        sigma_L = 0;
    else if (abs(minmod_a_L) < abs(minmod_b_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_a_L;
    else if (abs(minmod_b_L) < abs(minmod_a_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_b_L;

    return sigma - sigma_L;
}

//Van Albada 1 slope limiter
double slopeLimiter_vanAlbada1(double *arrayTemp, int j)
{
    double sigma;
    double sigma_L;
    double MC_ratio;
    double MC_ratio_L;

    MC_ratio = (arrayTemp[j] - arrayTemp[j-1])/(arrayTemp[j+1] - arrayTemp[j]);
    MC_ratio_L = (arrayTemp[j-1] - arrayTemp[j-2])/(arrayTemp[j] - arrayTemp[j-1]);
    
    sigma = (pow(MC_ratio,2) + MC_ratio)/(pow(MC_ratio,2) + 1);
    sigma_L = (pow(MC_ratio_L,2) + MC_ratio_L)/(pow(MC_ratio_L,2) + 1);

    return sigma - sigma_L;
}

//Lax Wendoff slope limiter
double slopeLimiter_LaxWendoff(double *arrayTemp, int j)
{
    double sigma;
    double sigma_L;
    return 0;
}