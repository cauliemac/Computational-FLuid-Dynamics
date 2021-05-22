#pragma once

/*
A header file containing slope limiters for 
use in differential equations
//TODO MOVE SLOPE LIMITERS TO A SEPERATE HEADER FILE
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

//size of grid
#define gridSize 100

//structs to hold the cell data
typedef struct CELL_STATE
{
    double Density[gridSize];
    double Pressure[gridSize];
    double Velocity[gridSize];

}cell_state;

typedef struct INTERFACE_CELL_STATE
{
    double Density;
    double Pressure;
    double Velocity;

}interface_cell_state;

//TODO UPDATE TO WORK WITH STRUCTS
//returns a grid with a different height between a and b
static double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine);

//a function to pick between slope limiters
double chooseSlopeLimiter(double *tempArray, int j, int SlopeType);

//a function to impliment the godunov Scheme
void GodunovScheme(cell_state temp_cell_state, cell_state* solution_cell_state, int j, double dx, double dt, interface_cell_state riemann_cell_state);

//starts the process of the next evolution of the godunov scheme
double AllEvolutions(cell_state solution_cell_state, cell_state temp_cell_state, int evolutions, double courant, double dx, interface_cell_state riemann_cell_state);

//gets the variable time step
double getDT(cell_state temp_cell_state, double courant, double dx);

// C function to find maximum in arr[] of size n
double largest(double arr[], int n);

//Calculates the resolved state given the left and right states.
void adiflux(cell_state temp_cell_state, int Left, int Right, double dx, double dt, interface_cell_state* riemann_cell_state);

//This function calculates the wave speed for a wave connecting states with pressures, pi, p ahead and behind respectively
double wave(double p,double pi, double gamma);