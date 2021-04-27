#pragma once

/*
A header file containing slope limiters for 
use in differential equations

*/
//monotonized Central Slope Limiter
double slopeLimiter_MC(double *arrayTemp, int j);

//Minmod slope Limiter
double slopeLimiter_minmod(double *arrayTemp, int j);

//Van Albada 1 slope limiter
double slopeLimiter_vanAlbada1(double *arrayTemp, int j);

//Lax Wendoff slope limiter
double slopeLimiter_LaxWendoff(double *arrayTemp, int j);