#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SlopeLimiters.h"

/*
1D Invisid Burgers Equation using finite volume method
With implimentation of slope limiters

TODO add slopelimiters to header files

*/

// declare starting variables

#define gridSize 500 //size of grid
const double gridSpacing = 2.0 / (gridSize);   //grid spacing ( also h)
const double evolutions = 1000;  //number of evolutions
const double timestepSize = 0.002;  //size of each timestep ( also k)
double courant = timestepSize/ gridSpacing; //Cournant number for printout


//creates 3 arrays to hold the grid values
double arraySolution[gridSize];
double arrayTemp[gridSize];
double initialConditions[gridSize];

//declareing the functions
static double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine); 
//double slopeLimiter_MC(double *arrayTemp, int j);
double AllEvolutions(double *arraySolution, int evolutions, double courant, double gridSpacing);
double BurgersEquation(double *arrayTemp, int j, int k);


//gives an initial conditions array with a square or sine wave (if sine == 1)
double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine)
{
    //create and open a file in write mode to store the initial conditions
    FILE *initial_file = NULL;
    initial_file = fopen("BurgersEquation_1D_results/BurgersEquationInitial.txt", "w");

    //populates array with one wavelength sine wave
    if (sine == 1)
        for (int i=0; i<gridSize; i++)
        {
            initialConditions[i] = sin(2*i/(gridSize/M_PI));
            fprintf(initial_file, "%f\n", initialConditions[i] );
        }

    else
        //populates the area between a and b as percentages of the grid with height 2
       for (int i=0; i<gridSize; i++)
        {                
            if (i/gridSize > a/100 && i/gridSize < b/100)   
                initialConditions[i] = 2;
            else
                initialConditions[i] = 1;
            
            fprintf(initial_file, "%f\n", initialConditions[i] );
        }

    //memcpy to copy initial conditions onto the solution array
    memcpy(arraySolution, initialConditions, gridSize * sizeof(double));

    return *initialConditions;
}

//Calls the Monotonized Central Slope Limiter
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

//k is to take the place of j+1 for upwind negative value (downwind) sections
double BurgersEquation(double *arrayTemp, int j, int k)
{
    double solution = arrayTemp[j] - courant * (0.5 * pow(arrayTemp[k],2) - 0.5*pow(arrayTemp[k-1],2)) - 0.5 * courant * (gridSpacing - timestepSize)*(slopeLimiter_MC(arrayTemp,k));
    return solution;
}

double AllEvolutions(double *arraySolution, int evolutions, double courant, double gridSpacing)
{
    for (int i = 0; i < evolutions; i++)
    {
        //copies the solutions array onto the temp array
        memcpy(arrayTemp, arraySolution, gridSize*sizeof(double));  

        //changes file name with evolution cycle.
        FILE *fpointer = NULL;
        char buffer[256]; // The filename buffer.

        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "BurgersEquation_1D_results/BurgersEquationSolution%i.txt", i);
        fpointer = fopen(buffer, "w");

        //calculates the next value of the current cell
        for (int j = 2; j < gridSize-1; j++)
        {
            //For the size of the grid, calculate the result of the 2nd order Burgers Equation
            if (arraySolution[j] > 0)

                // calculates next value of the array solution
                arraySolution[j] = BurgersEquation(arrayTemp, j, j);

            //for negative velocity
            //TODO this needs fixed. not calculating downwind correctly
            else
                //calculates next value of the array solution
                arraySolution[j] = BurgersEquation(arrayTemp, j, j+1);
            
            //uncomment below to print values to console
            //printf("%f\n", arraySolution[j])
            
            //print the x axis label (which is j) and the solution to a text file
            fprintf(fpointer, "%i \t %f\n", j, arraySolution[j]);
        }
    //printf("current evolutions is %i", i);
    }
}


int main ()
{   
    //calls the initial conditions function
    getInitialConditions(initialConditions, gridSize, 30, 60, 1);

    //call AllEvolutions to run
    AllEvolutions(arraySolution, evolutions, courant, gridSpacing);
   
    //print some useful info to the console
    printf("grid Spacing: %f\n", gridSpacing);
    printf("time step size: %f\n", timestepSize);
    printf("evolutions: %f\n", evolutions);
    printf("Courant number: %f\n", courant);

    return EXIT_SUCCESS;
}