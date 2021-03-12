/*
A header file containing slope limiters for 
use in differential equations

*/



void minmod(double a, double b, double c)
{   
    /*
    a is j-1
    b is j
    c is j+1
    */

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

    minmod_a_L = (arrayTemp[j-1] - arrayTemp[j-2])/gridSpacing;
    minmod_b_L = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing;

    if(minmod_a_L * minmod_b_L < 0)
        sigma_L = 0;
    else if (abs(minmod_a_L) < abs(minmod_b_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_a_L;
    else if (abs(minmod_b_L) < abs(minmod_a_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_b_L;

   /*
    double minmod_a;
    double minmod_b;
    double minmod_b_L;
    double minmod_a_L;
    double sigma;
    double sigma_L;

    minmod_a = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing;
    minmod_b = (arrayTemp[j+1] - arrayTemp[j])/gridSpacing;

    if(minmod_a * minmod_b < 0)
        sigma = 0;
    else if (abs(minmod_a) < abs(minmod_b) && minmod_a * minmod_b > 0)
        sigma = minmod_a;
    else if (abs(minmod_b) < abs(minmod_a) && minmod_a * minmod_b > 0)
        sigma = minmod_b;

    minmod_a_L = (arrayTemp[j-1] - arrayTemp[j-2])/gridSpacing;
    minmod_b_L = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing;

    if(minmod_a_L * minmod_b_L < 0)
        sigma_L = 0;
    else if (abs(minmod_a_L) < abs(minmod_b_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_a_L;
    else if (abs(minmod_b_L) < abs(minmod_a_L) && minmod_a_L * minmod_b_L > 0)
        sigma_L = minmod_b_L;
        */

}