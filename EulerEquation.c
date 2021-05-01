#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <Windows.h>
#include "functions.h"

/**************************************************************
 *  1D implimentation of Euler's Equations of Gas Dynamics    *
 *                                                            *
 *  g++ EulerEquation.c SlopeLimiters.c -o EulerEquation      *
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
float gamma_val = 5/3;
const int slope_limiter_type = 1   //1 for MC, 2 for Minmod, 3 for Van Albada 1

/*
//testing out structures for holding the data
typedef struct
{
    double solutionPressure[gridSize];
    double solutionDensity[gridSize];
    double solutionVelocity[gridSize];
}solution_cell_state;

typedef struct
{
    double tempPressure[gridSize];
    double tempDensity[gridSize];
    double tempVelocity[gridSize];
}temp_cell_state;
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

double chooseSlopeLimiter(double *tempArray, int j, int SlopeType);

double EulerEquationDensity(double *tempDensity, int j);
double EulerEquationMomentum(double *tempDensity, double *tempMomentum, double *tempEnergy, int j);
double EulerEquationEnergy(double *tempDensity, double *tempMomentum, double *tempEnergy, int j);
double pressure(double *tempDensity, double *tempMomentum, double *tempEnergy, int j);

double RiemannSolver(double *tempDensity, double *tempMomentum, double *tempEnergy, int Scheme, int j, int k);

double GodunovScheme(double *tempDensity, double *tempMomentum, double *tempEnergy, int j, int Scheme);

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
   
    for(int i=0; i<=100; ++i)
    {
        printf("\r[%3d%%]",i);
    } 
    printf("\n");
    system("pause");
    return 0;
}

/*
 *Returns an array with the values of one sine wave (if sine==1)
 *or a square wave with the peak between a and b as percentages
 *of the whole grid (if sine==0)
 */
//TODO return error if initial conditions are wrong
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
                initialConditions[i] = 0.5;
            }
            else
            {
                initialConditions[i] = 0.2;
            }
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            //fprintf(initial_momentum, " %i \t %f\n", i, initialConditions[i]);
            //fprintf(initial_energy, " %i \t %f\n", i, initialConditions[i]);
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
//TODO seperate opening files to a different funtion
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
            solutionDensity[j] = GodunovScheme(solutionDensity, solutionMomentum, solutionEnergy, j, 1);
            solutionMomentum[j] = GodunovScheme(solutionDensity, solutionMomentum, solutionEnergy, j, 2);
            solutionEnergy[j] = GodunovScheme(solutionDensity, solutionMomentum, solutionEnergy, j, 3);

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
double GodunovScheme(double *tempDensity, double *tempMomentum, double *tempEnergy, int j, int Scheme_DME)
{
    double Godunov;
    double densityLeft; double densityRight;
    double momentumLeft; double momentumRight;
    double energyLeft; double energyRight;



    if (Scheme_DME == 1)
    {
        densityLeft = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j-1, j);
        densityRight = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j, j+1);
        adiflux(double *tempDensity, double *tempPressure, double *tempVelocity, int Left, int Right, double dx, double dt, double *resolved_state)

        //densityLeft = adiflux(j-1, j, gridSpacing, timestepSize, j, 0, 0)
        
        if (tempDensity[j] > 0)
        {
            Godunov = tempDensity[j] - courant * (densityRight - densityLeft) - 0.5 * courant * (gridSpacing - timestepSize)*(chooseSlopeLimiter(tempDensity,j,slope_limiter_type));
        }
        else
        {
            Godunov = 1.0e-06;
            printf("Error in Riemann Density calculations at j=%i\n",j);
        }
        
    }

    else if (Scheme_DME == 2)
    {
        momentumLeft = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j-1, j);
        momentumRight = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j, j+1);

        if (tempMomentum[j] > 0)
        {
            Godunov = tempMomentum[j] - courant * (momentumRight - momentumLeft) - 0.5 * courant * (gridSpacing - timestepSize)*(chooseSlopeLimiter(tempMomentum,j,slope_limiter_type));
        }
        else
        {
            Godunov = 1.0e-06;
            printf("Error in Riemann Momentum calculations at j=%i\n",j);
        }
        
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
        energyLeft = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j-1, j);
        energyRight = RiemannSolver(tempDensity, tempMomentum, tempEnergy, Scheme_DME, j, j+1);

        if (tempEnergy[j] > 0)
        {
            Godunov = tempEnergy[j] - courant * (energyRight - energyLeft) - 0.5 * courant * (gridSpacing - timestepSize)*(chooseSlopeLimiter(tempEnergy,j,slope_limiter_type));
        }
        else
        {
            Godunov = 1.0e-06;
            printf("Error in Riemann Energy calculations at j=%i\n",j);
            //negative energy is not a physical solution (this case still happens cos of my shit code tho)
        }
        //printf("Godunov at Energy at %i is = %f\n",j,Godunov);
        //Sleep(500);
    }
    
    return Godunov;
}

//picks a slope limiter from a list in SlopeLimiters.c
double chooseSlopeLimiter(double *tempArray, int j, int SlopeType)
{
    double sigma;
    if (SlopeType == 1)
    {
        sigma = slopeLimiter_MC(tempArray,j);
    }
    else if(SlopeType == 2)
    {
        sigma = slopeLimiter_minmod(tempArray,j);
    }
    else if(SlopeType == 3)
    {
        sigma = slopeLimiter_vanAlbada1(tempArray,j);
    }
    return sigma;
}

/*
 *Returns the exact solution to the Riemann problem between two piecewise states
 *for the Euler Equations
 *chooses between the equations for Denstiy, Momentum, Energy
 *Need to use exact_adiabatic.c to find the solution for Eulers
 */
//TODO use adiflux() function from exact_adiabatic.c here for Riemann Solver 
double RiemannSolver(double *tempDensity, double *tempMomentum, double *tempEnergy, int Scheme, int Left, int Right)
{
    double Riemann;
    double EulerLeft; double EulerRight;

    /*
     *Changes the variables for the left and right state depending on 
     *whether it's for density, momentum, or energy
     */
    if (Scheme == 1)
    {
        EulerLeft = EulerEquationDensity(tempDensity, Left);
        EulerRight = EulerEquationDensity(tempDensity, Right);

        if (tempDensity[Left] >= tempDensity[Right])
        {
            Riemann = fmax(EulerLeft, EulerRight);
        }
        else if (tempDensity[Left] <= 0 && tempDensity[Right] >= 0)
        {
            Riemann = 0;
        }
        else
        {
            Riemann = fmin(EulerLeft, EulerRight);
        }

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    else if (Scheme == 2)
    {
        EulerLeft = EulerEquationMomentum(tempDensity, tempMomentum, tempEnergy, Left);
        EulerRight = EulerEquationMomentum(tempDensity, tempMomentum, tempEnergy, Right);

        if (tempMomentum[Left] >= tempMomentum[Right])
        {
            Riemann = fmax(EulerLeft, EulerRight);
        }
        else if (tempMomentum[Left] <= 0 && tempMomentum[Right] >= 0)
        {
            Riemann = 0;
        }
        else
        {
            Riemann = fmin(EulerLeft, EulerRight);
        }
        

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    else if (Scheme == 3)
    {
        EulerLeft = EulerEquationEnergy(tempDensity, tempMomentum, tempEnergy, Left);
        EulerRight = EulerEquationEnergy(tempDensity, tempMomentum, tempEnergy, Right);

        if (tempEnergy[Left] >= tempEnergy[Right])
        {
            Riemann = fmax(EulerLeft, EulerRight);
        }
        else if (tempEnergy[Left] <= 0 && tempEnergy[Right] >= 0)
        {
            Riemann = 0;
        }
        else
        {
            Riemann = fmin(EulerLeft, EulerRight);
        }

        //memccpy(array_for_Riemann, arrayTemp, gridSize*sizeof(double));
    }
    



    return Riemann;
}

/*
 *Calculates the solution to the Euler Equation at a given gridpoint
 *for the density
 *This is used to calculate the flux though a cell wall
 */
double EulerEquationDensity(double *tempDensity, int j)
/*
 * rho * u
 * which is q(2) in leveque, or tempDensity
 */
{
    double solution = tempDensity[j];

    return solution;
}

double EulerEquationMomentum(double *tempDensity, double *tempMomentum, double *tempEnergy, int j)
/*
 * (rho *  u^2) + pressure
 * which is q(2)^2 / q(1) + pressure(q)
 */
{
    double solution = (pow(tempMomentum[j],2) / tempDensity[j]) + pressure(tempEnergy, tempMomentum, tempDensity, j);

    return solution;
}

double EulerEquationEnergy(double *tempDensity, double *tempMomentum, double *tempEnergy, int j)
/*
 * (E + pressure) * u
 *
 * E = (pressure/gamma-1) + 1/2 * rho * u^2
 * 
 * which is q(2) * [q(3) pressure(q)] / q(1)
 */
{
    double solution = tempMomentum[j] * (tempEnergy[j] + pressure(tempDensity, tempMomentum, tempEnergy, j)) / tempDensity[j];

    return solution;
}

double pressure(double *tempDensity, double *tempMomentum, double *tempEnergy, int j)
{
    double pressure = tempEnergy[j] * (gamma_val-1) - 0.5 * (pow(tempMomentum[j],2) / tempDensity[j]);

    return pressure;
}