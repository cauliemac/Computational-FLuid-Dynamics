#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <Windows.h>
#include "SlopeLimiters.h"

/*
1D Invisid Burgers Equation using finite volume method
With implimentation of slope limiters

g++ BurgersEquation.c SlopeLimiters.c -o BurgersEquation
use line above to compile from multiple .c files
*/

// declare starting variables

#define gridSize 100 //size of grid
const double gridSpacing = 2.0 / (gridSize);   //grid spacing ( also h)
const double evolutions = 100;  //number of evolutions
const double timestepSize = 0.005;  //size of each timestep ( also k)
double courant = timestepSize/gridSpacing; //Cournant number for printout
//not exactly courant number, should be (max wave speed * timestep)/ gridSpacing


//creates 3 arrays to hold the grid values
double arraySolution[gridSize];
double arrayTemp[gridSize];
double initialConditions[gridSize];

//declareing the functions
static double getInitialConditions(double *initialConditions, int grid, float a, float b, int sine); 
double AllEvolutions(double *arraySolution, int evolutions, double courant, double gridSpacing);
double BurgersEquation(double *arrayTemp, int j);
double RiemannSolver(double *arrayTemp, int j, int k);
double GodunovScheme(double *arrayTemp, int j);

//TODO THIS DOES NOT WORK FOR SINE == 0, I.E. IF I DONT WANT A SINE WAVE
//gives an initial conditions array with a square or sine wave (if sine == 1)
double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine)
{
    //create and open a file in write mode to store the initial conditions
    FILE *initial_file = NULL;
    initial_file = fopen("BurgersEquation_1D_results/BurgersEquationInitial.txt", "w");

    //populates array with one wavelength sine wave
    if (sine == 1)
    {
        printf("Initial Conditions = Sine Wave");
        for (int i=0; i<gridSize; i++)
        {
            initialConditions[i] = sin(2*i/(gridSize/M_PI));
            fprintf(initial_file, " %i \t %f\n", i, initialConditions[i]);
        }
    }
    else
    {
        //populates the area between a and b as percentages of the grid with height 2
        printf("Initial Conditions = Square Wave from %d to %d", a, b);
        for (int i=0; i<gridSize; i++)
        {                
            float a_ratio = 100*a/gridSize;
            float b_ratio = 100*b/gridSize;

            if (i >= a_ratio && i <= b_ratio)   
            {
                initialConditions[i] = 2;
            }
            else
            {
                initialConditions[i] = 1;
            }
            fprintf(initial_file, " %i \t %f\n", i, initialConditions[i]);
        }
    }
    //memcpy to copy initial conditions onto the solution array
    memcpy(arraySolution, initialConditions, gridSize * sizeof(double));

    return *initialConditions;
}


//Take if pos/neg statement from all evolutions
double RiemannSolver(double *arrayTemp, int Left, int Right)
{
    double Riemann;
    double BurgerLeft = BurgersEquation(arrayTemp, Left);
    double BurgerRight = BurgersEquation(arrayTemp, Right);

    if (arrayTemp[Left] >= arrayTemp[Right])
    {
        Riemann = fmax(BurgerLeft, BurgerRight);
    }
    
    else if (arrayTemp[Left] <= 0 && arrayTemp[Right] >= 0)
    {
        Riemann = 0;
    }
    
    else
    {
        Riemann = fmin(BurgerLeft, BurgerRight);
    }

    return Riemann;
}

//k is to take the place of j+1 for upwind negative value (downwind) sections
//for upwind it is equal to j
double BurgersEquation(double *arrayTemp, int j)
{
    double solution = 0.5*pow(arrayTemp[j],2);

    return solution;
}

double GodunovScheme(double *arrayTemp, int j)
{
    double LeftBoundary = RiemannSolver(arrayTemp, j-1, j);
    double RightBoundary = RiemannSolver(arrayTemp, j, j+1);

    double Godunov = arrayTemp[j] - courant * (RightBoundary - LeftBoundary);//- 0.5 * courant * (gridSpacing - timestepSize)*(slopeLimiter_MC(arrayTemp,j));

    return Godunov;
}

//picks a slope limiter from a list in SlopeLimiters.c
//TODO get slope limiter picker working
/*
int chooseSlopeLimiter(int n);
{
    if (n = 1)


}
*/

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
            arraySolution[j] = GodunovScheme(arraySolution,j);

            //print the x axis label (which is j) and the solution to a text file
            fprintf(fpointer, "%i \t %f\n", j, arraySolution[j]);
            
        }
        fclose(fpointer);
    }
    return 0;
}

int main ()
{   
    //print some useful info to the console
    printf("grid Spacing: %f\n", gridSpacing);
    printf("time step size: %f\n", timestepSize);
    printf("evolutions: %f\n", evolutions);
    printf("Courant number: %f\n", courant);
    
    //calls the initial conditions function
    getInitialConditions(initialConditions, gridSize, 30, 60, 1);

    Sleep(2000);

    //call AllEvolutions to run
    AllEvolutions(arraySolution, evolutions, courant, gridSpacing);
   
    return 0;
}