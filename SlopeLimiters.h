#ifndef SLOPELIMITERS_H_   /* Include guard */
#define SLOPELIMITERS_H_

/*
A header file containing slope limiters for 
use in differential equations

*/

double slopeLimiter_MC(double *arrayTemp, int j);
/*
void minmod(double a, double b, double c, double d)
{   
    /*
    a is j-1
    b is j
    c is j+1
    d is j-2
    

    double minmod_a;
    double minmod_b;
    double minmod_b_L;
    double minmod_a_L;
    double sigma;
    double sigma_L;
    extern double gridSpacing;

    minmod_a = (b - a)/gridSpacing;
    minmod_b = (c - b)/gridSpacing;

    if(minmod_a * minmod_b < 0)
        sigma = 0;
    else if (abs(minmod_a) < abs(minmod_b) && minmod_a * minmod_b > 0)
        sigma = minmod_a;
    else if (abs(minmod_b) < abs(minmod_a) && minmod_a * minmod_b > 0)
        sigma = minmod_b;

    minmod_a_L = (a - arrayTemp[j-2])/gridSpacing;
    minmod_b_L = (b - a)/gridSpacing;

    if(minmod_a_L * minmod_b_L < 0)
        sigma_L = 0;
    else if (abs(minmod_a_L) < abs(minmod_b_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_a_L;
    else if (abs(minmod_b_L) < abs(minmod_a_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_b_L;

    
}
*/
#endif