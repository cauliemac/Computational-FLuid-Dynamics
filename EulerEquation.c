#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <Windows.h>
#include "SlopeLimiters.h"

/**************************************************************
 *  1D implimentation of Euler's Equations of Gas Dynamics    *
 *                                                            *
 *  g++ BurgersEquation.c SlopeLimiters.c -o BurgersEquation  *
 *  use line above to compile from multiple .c files          *
 *                                                            *
 **************************************************************/

// declare starting variables
#define gridSize 1000 //size of grid
const double gridSpacing = 2.0 / (gridSize);   //grid spacing ( also h)
const int evolutions = 100;  //number of evolutions
const float timestepSize = 0.005;  //size of each timestep ( also k)
double courant = timestepSize/gridSpacing; //Courant number for printout
/*
 *double courant is not exactly courant number, 
 *should be (wave speed * timestep)/ gridSpacing
 */


//creates 3 solution arrays to hold the grid values for Mass, Momentum, and Energy
double solutionMass[gridSize];
double solutionMomentum[gridSize];
double solutionEnergy[gridSize];

//creates three temp holding arrays for Mass, Momentum, and Energy
double tempMass[gridSize];
double tempMomentum[gridSize];
double tempEnergy[gridSize];

//creates three initial conditions arrays for Mass, Momentum, and Energy
double initialMass[gridSize];
double initialMomentum[gridSize];
double initialEnergy[gridSize];

//declaring the functions
static double getInitialConditions(double *initialConditions, int grid, float a, float b, int sine); 
double AllEvolutions(double *arraySolution, int evolutions, double courant, double gridSpacing);
double EulerEquation(double *arrayTemp, int j);
double RiemannSolver(double *arrayTemp, int j, int k);
double GodunovScheme(double *arrayTemp, int j);

/*
 *Returns an array with the values of one sine wave (if sine==1)
 *or a square wave with the peak between a and b as percentages
 *of the whole grid (if sine==0)
 */
double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine)
{
    //create and open a file in write mode to store the initial conditions
    FILE *initial_file = NULL;
    initial_file = fopen("EulerEquation_1D_results/EulerEquationInitial.txt", "w");

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
    //populates the area between a and b (as percentages of the grid) with height 2
    //TODO In future versions this height should be a variable
    else
    {
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
    fclose(initial_file);

    return *initialConditions;
}


/*
 *Returns the exact solution to the Riemann problem between two piecewise states
 *Currently set for Burgers equation
 *Need to use exact_adiabatic.c to find the solution for Eulers
 */
double RiemannSolver(double *arrayTemp, int Left, int Right)
{
    double Riemann;
    double EulerLeft = EulerEquation(arrayTemp, Left);
    double EulerRight = EulerEquation(arrayTemp, Right);

    if (arrayTemp[Left] >= arrayTemp[Right])
    {
        Riemann = fmax(EulerLeft, EulerRight);
    }
    
    else if (arrayTemp[Left] <= 0 && arrayTemp[Right] >= 0)
    {
        Riemann = 0;
    }
    
    else
    {
        Riemann = fmin(EulerLeft, EulerRight);
    }

    return Riemann;
}

/*
 *Calculates the solution to the Euler Equation at a given gridpoint
 *This is used to calculate the flux though a cell wall
 */
double EulerEquation(double *arrayTemp, int j)
{
    double solution = 0.5*pow(arrayTemp[j],2);

    return solution;
}

/*
 *Uses the Godunov scheme in combination with the flux in and out of a gridpoint
 *to calculate the next value of a gridpoint
 */
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

/*
 *Copies the solution array to the temp array
 *calls the Godunov Scheme function for all evolutions
 *and prints the values to a text file
 *closes text file
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
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerEquationSolution%i.txt", i);
        fpointer = fopen(buffer, "w");

        //calculates the next value of the current cell
        for (int j = 2; j < gridSize-1; j++)
        {
            arraySolution[j] = GodunovScheme(arraySolution,j);
            float x = arraySolution[j];
            printf("Array Solution at %i is: %f\n", j, x);

            //print the x axis label (which is j) and the solution to a text file
            fprintf(fpointer, "%i \t %f\n", j, arraySolution[j]);
            
        }
        fclose(fpointer);
    }
    return 0;
}

/*
 *Prints some useful info to the console
 *Gets the initial conditions
 *Calls the AllEvolutions function
 */
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