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

//Newton-Raphson support function for calculating f and df
void func_newt(double x,double *f, double *df, double t1, double t2, double t3, double *gl, double *gr);

//Newton-Raphson iterator
double newt(double x, double xmin, double xmax, double t1, double t2, double t3, double *gl, double *gr, int *newt_check);

//Function evaluation for bisection method iterator
double func_bis(double x, double t1, double t2, double t3, double *gl, double *gr);

//Bisection method iterator
double rtbis(double x1,double x2,double xacc, double t1, double t2, double t3, double *gl, double *gr);

//Calculates various functions of gamma, the ratio of specific heats
void gamma_calc(double gamma, double *g);

//Calculates the resolved state given the left and right states. 
void adiflux(solution_cell_state, double dx, double dt, double *resolved_state);