#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <Windows.h>
#include "SlopeLimiters.h"

/**************************************************************
 *  1D implimentation of Euler's Equations of Gas Dynamics    *
 *                                                            *
 *  g++ EulerEquation.c SlopeLimiters.c -o EulerEquation  *
 *  use line above to compile from multiple .c files          *
 *                                                            *
 **************************************************************/

// declare starting variables
#define gridSize 100 //size of grid
const double gridSpacing = 2.0 / (gridSize);   //grid spacing ( also h)
const int evolutions = 100;  //number of evolutions
const float timestepSize = 0.005;  //size of each timestep ( also k)
double courant = timestepSize/gridSpacing; //Courant number for printout
/*
 *double courant is not exactly courant number, 
 *should be (wave speed * timestep)/ gridSpacing
 */

//creates 3 solution arrays to hold the grid values for Mass, Momentum, and Energy
double solutionDensity[gridSize];
double solutionMomentum[gridSize];
double solutionEnergy[gridSize];

//creates three temp holding arrays for Mass, Momentum, and Energy
double tempDensity[gridSize];
double tempMomentum[gridSize];
double tempEnergy[gridSize];

//creates three initial conditions arrays for Mass, Momentum, and Energy
double initialDensity[gridSize];
double initialMomentum[gridSize];
double initialEnergy[gridSize];

//declaring the functions
double initialConditions[gridSize];
static double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine); 
double AllEvolutions(double *solutionDensity, double *solutionMomentum, double *solutionEnergy, int evolutions, double courant, double gridSpacing);

double EulerEquationDensity(double *tempDensity, int j);
double EulerEquationMomentum(double *tempMomentum, int j);
double EulerEquationEnergy(double *tempEnergy, int j);

double RiemannSolver(double *arrayTemp, int Scheme, int j, int k);

double GodunovScheme(double *arrayTemp, int j, int Scheme);

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
    getInitialConditions(initialConditions, gridSize, 0, 50, 0);

    Sleep(2000);

    //call AllEvolutions to run
    AllEvolutions(solutionDensity, solutionMomentum, solutionEnergy, evolutions, courant, gridSpacing);
   
    return 0;
}

/*
 *Returns an array with the values of one sine wave (if sine==1)
 *or a square wave with the peak between a and b as percentages
 *of the whole grid (if sine==0)
 */
static double getInitialConditions(double *initialConditions, int grid, int a, int b, int sine)
{
    //create and open a file in write mode to store the initial conditions
    FILE *initial_density = NULL;
    initial_density = fopen("EulerEquation_1D_results/_EulerInitialDensity.txt", "w");
    FILE *initial_momentum = NULL;
    initial_momentum = fopen("EulerEquation_1D_results/_EulerInitialMomentum.txt", "w");
    FILE *initial_energy = NULL;
    initial_energy = fopen("EulerEquation_1D_results/_EulerInitialEnergy.txt", "w");

    //populates array with one wavelength sine wave
    if (sine == 1)
    {
        printf("Initial Conditions = Sine Wave");
        for (int i=0; i<gridSize; i++)
        {
            initialConditions[i] = sin(2*i/(gridSize/M_PI));
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_momentum, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_energy, " %i \t %f\n", i, initialConditions[i]);
        }
    }
    //populates the area between a and b (as percentages of the grid) with height 2
    //TODO In future versions this height should be a variable
    else
    {
        printf("Initial Conditions = Square Wave from %d to %d", a, b);
        for (int i=0; i<gridSize; i++)
        {                
            float a_ratio = a*gridSize/100;
            float b_ratio = b*gridSize/100;

            if (i >= a_ratio && i <= b_ratio)   
            {
                initialConditions[i] = 2;
            }
            else
            {
                initialConditions[i] = 1;
            }
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_momentum, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_energy, " %i \t %f\n", i, initialConditions[i]);
        }
    }
    //memcpy to copy initial conditions onto the solution array
    memcpy(solutionDensity, initialConditions, gridSize * sizeof(double));
    memcpy(solutionMomentum, initialConditions, gridSize * sizeof(double));
    memcpy(solutionEnergy, initialConditions, gridSize * sizeof(double));
    
    fclose(initial_density);
    fclose(initial_momentum);
    fclose(initial_energy);

    return 0;
}

/*
 *Copies the solution array to the temp array
 *calls the Godunov Scheme function for all evolutions
 *and prints the values to a text file
 *closes text file
 */
double AllEvolutions(double *solutionDensity, double *solutionMomentum, double *solutionEnergy, int evolutions, double courant, double gridSpacing)
{
    for (int i = 0; i < evolutions; i++)
    {
        //copies the solutions array onto the temp array
        memcpy(tempDensity, solutionDensity, gridSize*sizeof(double));
        memcpy(tempMomentum, solutionMomentum, gridSize*sizeof(double));
        memcpy(tempEnergy, solutionEnergy, gridSize*sizeof(double)); 

        /*
         *Creates a text file for density, momentum, and energy wwith each evolution
         *changes file name with evolution cycle.
         */
        FILE *densityFile = NULL;
        char buffer[256]; // The filename buffer.
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerDensitySolution%i.txt", i);
        densityFile = fopen(buffer, "w");

        FILE *momentumFile = NULL;
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerMomentumSolution%i.txt", i);
        momentumFile = fopen(buffer, "w");

        FILE *energyFile = NULL;
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerEnergySolution%i.txt", i);
        energyFile = fopen(buffer, "w");

        //calculates the next value of the current cell
        for (int j = 1; j < gridSize-1; j++)
        {
            /*
             *Calls the Godunov Scheme for Density, Momentum, and Energy
             *if Scheme == 1 then it uses the varialbes for density
             *if Scheme == 2 then it uses momentum and so on
             */
            solutionDensity[j] = GodunovScheme(solutionDensity, j, 1);
            solutionMomentum[j] = GodunovScheme(solutionMomentum, j, 2);
            solutionEnergy[j] = GodunovScheme(solutionEnergy, j, 3);

            //printf("Momentum at %i is =  %f\n", j, solutionMomentum[j]);
 
            //print the x axis label (which is j) and the solution to a text file
            fprintf(densityFile, "%i \t %f\n", j, solutionDensity[j]);
            fprintf(momentumFile, "%i \t %f\n", j, solutionMomentum[j]);
            fprintf(energyFile, "%i \t %f\n", j, solutionEnergy[j]);
            
        }
        fclose(densityFile);
        fclose(momentumFile);
        fclose(energyFile);
    }
    return 0;
}

/*
 *Uses the Godunov scheme in combination with the flux in and out of a gridpoint
 *to calculate the next value of a gridpoint

 *if Scheme == 1 then it uses the varialbes for density
 *if Scheme == 2 then it uses the variables for momentum
 *if Scheme == 3 the it uses the variables for Energy
 */
double GodunovScheme(double *arrayTemp, int j, int Scheme_DME)
{
    double Godunov;
    double densityLeft; double densityRight;
    double momentumLeft; double momentumRight;
    double energyLeft; double energyRight;



    if (Scheme_DME == 1)
    {
        densityLeft = RiemannSolver(arrayTemp, Scheme_DME, j-1, j);
        densityRight = RiemannSolver(arrayTemp, Scheme_DME, j, j+1);
        

        Godunov = arrayTemp[j] - courant * (densityRight - densityLeft);//- 0.5 * courant * (gridSpacing - timestepSize)*(slopeLimiter_MC(tempDensity,j));
    }

    else if (Scheme_DME == 2)
    {
        momentumLeft = RiemannSolver(arrayTemp, Scheme_DME, j-1, j);
        momentumRight = RiemannSolver(arrayTemp, Scheme_DME, j, j+1);

        Godunov = arrayTemp[j] - courant * (momentumRight - momentumLeft);
        /*
        printf("here\n");
        printf("value for arrayTemp at %i is = %f\n",j,arrayTemp[j]);
        printf("value for momentum is = %f\n",Godunov);
        printf("value for momentumLEFT is = %f\n", momentumLeft);
        printf("value for momentumRIGHT is = %f\n", momentumRight);
        printf("temp value for momentum is = %f\n", tempMomentum);
        printf("value for courant is = %f\n", courant);
        Sleep(1000);
        */
        
    }

    else if (Scheme_DME == 3)
    {
        energyLeft = RiemannSolver(arrayTemp, Scheme_DME, j-1, j);
        energyRight = RiemannSolver(arrayTemp, Scheme_DME, j, j+1);

        Godunov = arrayTemp[j] - courant * (energyRight - energyLeft);
        //printf("Godunov at Energy at %i is = %f\n",j,Godunov);
        //Sleep(500);
    }
    
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
 *Returns the exact solution to the Riemann problem between two piecewise states
 *for the Euler Equations
 *chooses between the equations for Denstiy, Momentum, Energy
 *Need to use exact_adiabatic.c to find the solution for Eulers
 */
//TODO use adiflux() function from exact_adiabatic.c here for Riemann Solver 
double RiemannSolver(double *arrayTemp, int Scheme, int Left, int Right)
{
    double Riemann;
    double EulerLeft; double EulerRight;

    /*
     *Changes the variables for the left and right state depending on 
     *whether it's for density, momentum, or energy
     */
    if (Scheme == 1)
    {
        EulerLeft = EulerEquationDensity(arrayTemp, Left);
        EulerRight = EulerEquationDensity(arrayTemp, Right);

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    else if (Scheme == 2)
    {
        EulerLeft = EulerEquationMomentum(arrayTemp, Left);
        EulerRight = EulerEquationMomentum(arrayTemp, Right);

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    else if (Scheme == 3)
    {
        EulerLeft = EulerEquationEnergy(arrayTemp, Left);
        EulerRight = EulerEquationEnergy(arrayTemp, Right);

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    

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
 *for the density
 *This is used to calculate the flux though a cell wall
 */
//TODO  Use Euler's Equations here
double EulerEquationDensity(double *tempDensity, int j)
{
    double solution = 0.5*pow(tempDensity[j],2);

    return solution;
}

double EulerEquationMomentum(double *tempMomentum, int j)
{
    double solution = 0.5*pow(tempMomentum[j],2);

    return solution;
}

double EulerEquationEnergy(double *tempEnergy, int j)
{
    double solution = 0.5*pow(tempEnergy[j],2);

    return solution;
}