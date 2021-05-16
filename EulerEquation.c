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
//#define gridSize 100 //size of grid
const double dx = 2.0 / (gridSize);   //grid spacing ( also h)
const int evolutions = 100;  //number of evolutions
const float dt = 0.005;  //size of each timestep ( also k)
double courant = dt/dx; //Courant number for printout
/*
 *double courant is not exactly courant number, 
 *should be (wave speed * timestep)/ dx
 */
float gamma_val = 5/3;
const int slope_limiter_type = 1;   //1 for MC, 2 for Minmod, 3 for Van Albada 1

cell_state temp_cell_state, solution_cell_state;
interface_cell_state riemann_cell_state; 

//creates three initial conditions arrays for Mass, Momentum, and Energy
double initialDensity[gridSize];
double initialPressure[gridSize];
double initialVelocity[gridSize];

//declaring the functions
double initialConditions[gridSize];

//Calculates the resolved state given the left and right states.
//void adiflux cell_state temp_cell_state, int Left, int Right, double dx, double dt, interface_cell_state riemann_cell_state);
//TODO redo this as a void 

/*
 *Prints some useful info to the console
 *Gets the initial conditions
 *Calls the AllEvolutions function
 */
int main ()
{   
    //print some useful info to the console
    printf("grid Spacing: %f\n", dx);
    printf("time step size: %f\n", dt);
    printf("evolutions: %f\n", evolutions);
    printf("Courant number: %f\n", courant);

    //calls the initial conditions function
    getInitialConditions(initialConditions, gridSize, 0, 15, 0);

    Sleep(2000);

    //call AllEvolutions to run
    AllEvolutions(solution_cell_state, temp_cell_state, evolutions, courant, dx, riemann_cell_state);
   
    for(int i=0; i<=100; ++i)
    {
        printf("\r[%3d%%]",i);
        //Sleep(10);
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
    FILE *initial_pressure = NULL;
    initial_pressure = fopen("EulerEquation_1D_results/_EulerInitialPressure.txt", "w");
    FILE *initial_velocity = NULL;
    initial_velocity = fopen("EulerEquation_1D_results/_EulerInitialVelocity.txt", "w");

    //populates array with one wavelength sine wave
    if (sine == 1)
    {
        printf("Initial Conditions = Sine Wave");
        for (int i=0; i<gridSize; i++)
        {
            initialConditions[i] = sin(2*i/(gridSize/M_PI));
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_pressure, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_velocity, " %i \t %f\n", i, initialConditions[i]);
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
            
            //writing files to solution text file
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_pressure, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_velocity, " %i \t %f\n", i, 0); //velocity is 0
        }
    }
    //memcpy to copy initial conditions onto the solution array
    memcpy(solution_cell_state.Density, initialConditions, gridSize * sizeof(double));
    memcpy(solution_cell_state.Pressure, initialConditions, gridSize * sizeof(double));
    memcpy(solution_cell_state.Velocity, initialConditions, gridSize * sizeof(double));
    
    fclose(initial_density);
    fclose(initial_pressure);
    fclose(initial_velocity);

    return 0;
}

/*
 *Copies the solution array to the temp array
 *calls the Godunov Scheme function for all evolutions
 *and prints the values to a text file
 *closes text file
 */
//TODO seperate opening files to a different funtion
double AllEvolutions(cell_state solution_cell_state, cell_state temp_cell_state, int evolutions, double courant, double dx, interface_cell_state riemann_cell_state){

    for (int i = 0; i < evolutions; i++)
    {
        //copies the solutions array onto the temp array
        temp_cell_state = solution_cell_state;
        //memcpy(temp_cell_state.Density, solution_cell_state.Density, gridSize * sizeof(double));
        //memcpy(temp_cell_state.Pressure, solution_cell_state.Pressure, gridSize * sizeof(double));
        //memcpy(temp_cell_state.Velocity, solution_cell_state.Velocity, gridSize * sizeof(double));

        /*
         *Creates a text file for density, momentum, and energy wwith each evolution
         *changes file name with evolution cycle.
         */
        FILE *densityFile = NULL;
        char buffer[256]; // The filename buffer.
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerDensitySolution%i.txt", i);
        densityFile = fopen(buffer, "w");

        FILE *pressureFile = NULL;
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerPressureSolution%i.txt", i);
        pressureFile = fopen(buffer, "w");

        FILE *velocityFile = NULL;
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerVelocitySolution%i.txt", i);
        velocityFile = fopen(buffer, "w");

        //calculates the next value of the current cell
        for (int j = 1; j < gridSize-1; j++)
        {
           
            GodunovScheme(temp_cell_state, &solution_cell_state, j, dx, dt, riemann_cell_state);

            //printf("Pressure at %i is =  %f\n", j, solutionVelocity[j]);
 
            //print the x axis label (which is j) and the solution to a text file
            //fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);

            //printf("Solution densisty after Godunov at j=%i is %f\n",j,solution_cell_state.Density[j]);
            //printf("-------------------------------------------------------\n");
            //system("pause");

            //solution_cell_state.Density[j] = 10;

            fprintf(densityFile, "%i \t %f\n", j, solution_cell_state.Density[j]);
            fprintf(pressureFile, "%i \t %f\n", j, solution_cell_state.Pressure[j]);
            fprintf(velocityFile, "%i \t %f\n", j, solution_cell_state.Velocity[j]);
        }
        /*for (int j = 1; j < gridSize-1; j++)
        {
           fprintf(densityFile, "%i test \t %f\n", j, solution_cell_state.Density[j]);
           fprintf(pressureFile, "%i \t %f\n", j, solution_cell_state.Pressure[j]);
           fprintf(velocityFile, "%i \t %f\n", j, solution_cell_state.Velocity[j]);
        }
        */
        fclose(densityFile);
        fclose(pressureFile);
        fclose(velocityFile);
    }
    return 0;
}

/*
 *Uses the Godunov scheme in combination with the flux in and out of a gridpoint
 *to calculate the next value of a gridpoint

 */
void GodunovScheme (cell_state temp_cell_state, cell_state* solution_cell_state, int j, double dx, double dt, interface_cell_state riemann_cell_state){

    double Godunov;
    double densityLeft, densityRight;
    double pressureLeft, pressureRight;
    double velocityLeft, velocityRight;
    int Left, Right;

    if (solution_cell_state->Velocity[j] >= 0)
    {
        Left = j-1;
        Right = j;

        /*
        *Returns the exact solution to the Riemann problem between 
        *two piecewise states for the Euler Equations.
        *chooses between the equations for pressure, density, and velocity
        *Need to use exact_adiabatic.c to find the solution for Eulers
        */
        adiflux(temp_cell_state, Left, Right, dx, dt, &riemann_cell_state);
        
        densityLeft = riemann_cell_state.Density;
        pressureLeft = riemann_cell_state.Pressure;
        velocityLeft = riemann_cell_state.Velocity;

        Left = j;
        Right = j+1;
        
        adiflux(temp_cell_state, Left, Right, dx, dt, &riemann_cell_state);

        /*
        printf("the adiflux return at j = %i is:\n",j);
        printf("Density is = %f\n",riemann_cell_state.Density);
        printf("Pressure is = %f\n",riemann_cell_state.Pressure);
        printf("Velocity is = %f\n",riemann_cell_state.Velocity);
        printf("\n");
        */

        //Sleep(2000);

        densityRight = riemann_cell_state.Density;
        pressureRight = riemann_cell_state.Pressure;
        velocityRight = riemann_cell_state.Velocity;
    }
    else
    {
        Left = j;
        Right = j+1;
        printf("\n ERROR IN NEGATIVE VELOCITY!!!!!\nEXITING");
        //TODO THIS NEEDS TO BE FIXED
        system("pause");
    }
 
    //TODO fix slope limiters for res

    //return solution_cell_state;
    solution_cell_state->Density[j] = temp_cell_state.Density[j] - courant * (densityRight - densityLeft);// - 0.5 * courant * (dx - dt)*(chooseSlopeLimiter(temp_cell_state.Density,j,slope_limiter_type));
    solution_cell_state->Pressure[j] = temp_cell_state.Pressure[j] - courant * (pressureRight - pressureLeft);
    solution_cell_state->Velocity[j] = temp_cell_state.Velocity[j] - courant * (velocityRight - velocityLeft);

    //printf("solution density at j=%i in Godunov is %f\n",j,solution_cell_state.Density[j]);
    //system("pause");

    double temppp;
    temppp = courant * (densityRight - densityLeft);

    
    /*
    printf("Density right at %d is = %f\n",j,densityRight);
    printf("Pressure right at %d is = %f\n",j,pressureRight);
    printf("Velocity right at %d is = %f\n",j,velocityRight);
    //system("\npause\n");
    printf("\n");

    printf("Density left at %d is = %f\n",j,densityLeft);
    printf("Pressure left at %d is = %f\n",j,pressureLeft);
    printf("Velocity left at %d is = %f\n",j,velocityLeft);
    system("\npause\n");
    printf("\n");
    /*
    if (j == 16)
    {
        printf("solution density is = %f   \n", solution_cell_state.Density[j]);
        printf("correction is = %f\n", temppp);
        
        printf("the adiflux return at j = %i is:\n",j);
        
        

        system("\npause\n");
    }
    */
}

//picks a slope limiter from a list in SlopeLimiters.c
//TODO impliment slopelimiter picker
/*
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
*/