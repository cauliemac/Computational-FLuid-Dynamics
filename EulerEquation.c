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
const double dx = 2.0 / (gridSize);   //grid spacing ( also h)
const int evolutions = 100;  //number of evolutions
double dt = 0.005;  //size of each timestep ( also k)
double courant = 0.2; //Desired Courant number for printout
/*
 *Cournant number is C = (wave speed * timestep)/ dx
 */
float gamma_val = 5/3;
const int slope_limiter_type = 1;   //1 for MC, 2 for Minmod, 3 for Van Albada 1
int file_skipper = 5; //only prints every n files

cell_state temp_cell_state, solution_cell_state;
interface_cell_state riemann_cell_state;
conservative_variables temp_conserve, solution_conserve;

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
    printf("The grid size is %f\n",gridSize);
    printf("The grid spacing is: %f\n", dx);
    printf("For %d evolutions\n", evolutions);
    printf("With a Courant number of %f\n", courant);

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
static double getInitialConditions(double *initialConditions, int grid, int a, int b, int is_sine)
{
    //create and open a file in write mode to store the initial conditions
    FILE *initial_density = NULL;
    initial_density = fopen("EulerEquation_1D_results/_EulerInitialDensity.txt", "w");
    FILE *initial_pressure = NULL;
    initial_pressure = fopen("EulerEquation_1D_results/_EulerInitialPressure.txt", "w");
    FILE *initial_velocity = NULL;
    initial_velocity = fopen("EulerEquation_1D_results/_EulerInitialVelocity.txt", "w");

    double init_velocity_mod[gridSize];
    
    //populates array with one wavelength sine wave
    if (is_sine == 1)
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
        printf("Initial Conditions = Square Wave from %d to %d percent of grid\n", a, b);
        for (int i=0; i<gridSize; i++)
        {                
            float a_ratio = a*gridSize/100;
            float b_ratio = b*gridSize/100;

            if (i >= a_ratio && i <= b_ratio)   
            {
                initialConditions[i] = 2.5;
            }
            else
            {
                initialConditions[i] = 2;
            }

            init_velocity_mod[i] = 0.5;
            
            //writing files to solution text file
            fprintf(initial_density, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_pressure, " %i \t %f\n", i, initialConditions[i]);
            fprintf(initial_velocity, " %i \t %f\n", i, init_velocity_mod[i]); //velocity is 0
        }
    }
    //memcpy to copy initial conditions onto the solution array
    memcpy(solution_cell_state.Density, initialConditions, gridSize * sizeof(double));
    memcpy(solution_cell_state.Pressure, initialConditions, gridSize * sizeof(double));
    memcpy(solution_cell_state.Velocity, init_velocity_mod, gridSize * sizeof(double));

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
        /*
         *Creates a text file for density, momentum, and energy wwith each evolution
         *changes file name with evolution cycle.
         *Puts "file" then, "i" then, ".txt" in to filename.
         *opens files in write mode
        */

        if(i%file_skipper == 0)
        {
            char buffer[256]; // The filename buffer.
            FILE *densityFile = NULL;
            snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerDensitySolution%i.txt", i);
            densityFile = fopen(buffer, "w");
            FILE *pressureFile = NULL;
            snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerPressureSolution%i.txt", i);
            pressureFile = fopen(buffer, "w");
            FILE *velocityFile = NULL; 
            snprintf(buffer, sizeof(char) * 256, "EulerEquation_1D_results/EulerVelocitySolution%i.txt", i);
            velocityFile = fopen(buffer, "w");
        }
        
        //copies the solutions array onto the temp array
        temp_cell_state = solution_cell_state;
        //temp_conserve = solution_conserve;
        
        //gets the variable time step dt
        dt = getDT(temp_cell_state, courant, dx);

        //sets the values for the conserved variables for the evolution
        for (int i = 0; i< gridSize; i++)
        {
            temp_conserve.Mass[i] = temp_cell_state.Density[i] * dx;
            temp_conserve.Momentum[i] = temp_cell_state.Density[i] * temp_cell_state.Velocity[i] * dx;
            temp_conserve.Energy[i] = ((temp_cell_state.Pressure[i]/(2.0/3.0)) + 0.5 * temp_cell_state.Density[i] * (pow(temp_cell_state.Velocity[i],2.0))) * dx;
        }

        //calculates the next value of the current cell
        for (int j = 1; j < gridSize-2; j++)
        {
            GodunovScheme(temp_cell_state, &solution_cell_state, temp_conserve, solution_conserve, j, dx, dt, riemann_cell_state);

            if(i%file_skipper == 0)
            {
                fprintf(densityFile, "%i \t %f\n", j, solution_cell_state.Density[j]);
                fprintf(pressureFile, "%i \t %f\n", j, solution_cell_state.Pressure[j]);
                fprintf(velocityFile, "%i \t %f\n", j, solution_cell_state.Velocity[j]);
            }

            
        }
        
        if(i%file_skipper == 0)
        {
            fclose(densityFile);
            fclose(pressureFile);
            fclose(velocityFile);
        }
    }
    return 0;
}

/*
 *Uses the Godunov scheme in combination with the flux in and out of a gridpoint
 *to calculate the next value of a gridpoint

 */
void GodunovScheme (cell_state temp_cell_state, cell_state* solution_cell_state, conservative_variables temp_conserve, conservative_variables solution_conserve, int j, double dx, double dt, interface_cell_state riemann_cell_state)
{

    double densityLeft, densityRight;
    double pressureLeft, pressureRight;
    double velocityLeft, velocityRight;
    int Left, Right;
    double sound_speed;
    double c;//signal speed
    double c_Left, c_right;

    double MassFluxLeft, MomentumFluxLeft, EnergyFluxLeft;
    double MassFluxRight, MomentumFluxRight, EnergyFluxRight;

    if (solution_cell_state->Velocity[j] >= 0)
    {
        //Left cell interfaces
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
        
        //convert primitive variables to conservative variables
        MassFluxLeft = densityLeft * velocityLeft * dx;
        MomentumFluxLeft = densityLeft * pow(velocityLeft,2.0) + pressureLeft * dx;
        EnergyFluxLeft = velocityLeft*(((3.0*pressureLeft /(2.0)) + 0.5 * densityLeft * (pow(velocityLeft,2.0))) + pressureLeft) * dx;

        //right cell interfaces
        Left = j;
        Right = j+1;
        
        adiflux(temp_cell_state, Left, Right, dx, dt, &riemann_cell_state);

        densityRight = riemann_cell_state.Density;
        pressureRight = riemann_cell_state.Pressure;
        velocityRight = riemann_cell_state.Velocity;
        
        //convert primitive variables to conservative variables
        MassFluxRight = densityRight * velocityRight * dx;
        MomentumFluxRight = densityRight * pow(velocityRight,2.0) + pressureRight * dx;
        EnergyFluxRight = velocityRight*(((pressureRight/(2.0/3.0)) + 0.5 * densityRight * (pow(velocityRight,2.0)))+ pressureRight) * dx;

        //calculate the sound speed and signal speed of the wave
        sound_speed = abs(sqrt(gamma_val*temp_cell_state.Pressure[j]/temp_cell_state.Density[j]));
        c = sound_speed * temp_cell_state.Velocity[j];
        //c = signal speed

        //Use Godunov method to update the value of the cell.
        solution_conserve.Mass[j] = temp_conserve.Mass[j] - (dt/dx) * (MassFluxRight - MassFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Mass,j,slope_limiter_type));
        solution_conserve.Momentum[j] = temp_conserve.Momentum[j] - (dt/dx) * (MomentumFluxRight - MomentumFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Momentum,j,slope_limiter_type));
        solution_conserve.Energy[j] = temp_conserve.Energy[j] - (dt/dx) * (EnergyFluxRight - EnergyFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Energy,j,slope_limiter_type));

        //convert back to primitave variables and return solution_cell_state
        solution_cell_state->Density[j] = solution_conserve.Mass[j] / dx;
        solution_cell_state->Velocity[j] = (solution_conserve.Momentum[j]/solution_cell_state->Density[j]) / dx;
        solution_cell_state->Pressure[j] = ((solution_conserve.Energy[j]/dx) - (0.5 * solution_cell_state->Density[j] * pow(solution_cell_state->Velocity[j],2.0))) * (2.0/3.0);
        
        /*
        printf("#######################################################\n");
        printf("for cell j=%d\n\n",j);
        
        printf("Density Left is %f\n",densityLeft);
        printf("Density Right is %f\n",densityRight);
        printf("velocity Left is %f\n",velocityLeft);
        printf("velocity Right is %f\n",velocityRight);
        printf("Pressure Left is %f\n",pressureLeft);
        printf("Pressure Right is %f\n\n",pressureRight);
        
        printf("Mass Left is %f\n",MassLeft);
        printf("Mass Right is %f\n",MassRight);
        printf("Momentum Left is %f\n",MomentumLeft);
        printf("Momentum Right is %f\n",MomentumRight);
        printf("Energy Left is %f\n",EnergyLeft);
        printf("Energy Right is %f\n\n",EnergyRight);

        printf("---------------------------\n");
        
        printf("solution Mass is %f\n",solution_conserve.Mass[j]);
        printf("solution Momentum is %f\n",solution_conserve.Momentum[j]);
        printf("solution energy is %f\n\n",solution_conserve.Energy[j]);
        
        printf("Density Solution is %f\n",solution_cell_state->Density[j]);
        printf("velocity Solution is %f\n",solution_cell_state->Velocity[j]);
        printf("Pressure Solution is %f\n\n",solution_cell_state->Pressure[j]);
        
        system("pause");
        */

        //Conservation of Momentum = rho*u
        //Conservation of Energy = (rho * u^2) + p
        //Equation of state = u*(E + p)
    }
    else
    {
        //Left cell interfaces
        Left = j;
        Right = j+1;

        adiflux(temp_cell_state, Left, Right, dx, dt, &riemann_cell_state);
        
        densityLeft = riemann_cell_state.Density;
        pressureLeft = riemann_cell_state.Pressure;
        velocityLeft = riemann_cell_state.Velocity;

        //convert primitive variables to conservative variables
        MassFluxLeft = densityLeft * velocityLeft * dx;
        MomentumFluxLeft = densityLeft * pow(velocityLeft,2.0) + pressureLeft * dx;
        EnergyFluxLeft = velocityLeft*(((3.0*pressureLeft /(2.0)) + 0.5 * densityLeft * (pow(velocityLeft,2.0))) + pressureLeft) * dx;

        //right cell interfaces
        Left = j+1;
        Right = j+2;
        
        adiflux(temp_cell_state, Left, Right, dx, dt, &riemann_cell_state);

        densityRight = riemann_cell_state.Density;
        pressureRight = riemann_cell_state.Pressure;
        velocityRight = riemann_cell_state.Velocity;

        //convert primitive variables to conservative variables
        MassFluxRight = densityRight * velocityRight * dx;
        MomentumFluxRight = densityRight * pow(velocityRight,2.0) + pressureRight * dx;
        EnergyFluxRight = velocityRight*(((pressureRight/(2.0/3.0)) + 0.5 * densityRight * (pow(velocityRight,2.0)))+ pressureRight) * dx;

        //calculate the sound speed and signal speed of the wave
        sound_speed = abs(sqrt(gamma_val*temp_cell_state.Pressure[j]/temp_cell_state.Density[j]));
        c = sound_speed * temp_cell_state.Velocity[j];
        //c = signal speed

        //Use Godunov method to update the value of the cell.
        solution_conserve.Mass[j] = temp_conserve.Mass[j] - (dt/dx) * (MassFluxRight - MassFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Mass,j,slope_limiter_type));
        solution_conserve.Momentum[j] = temp_conserve.Momentum[j] - (dt/dx) * (MomentumFluxRight - MomentumFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Momentum,j,slope_limiter_type));
        solution_conserve.Energy[j] = temp_conserve.Energy[j] - (dt/dx) * (EnergyFluxRight - EnergyFluxLeft);// - 0.5 * (c*dt/dx) * (dx - c*dt)*(chooseSlopeLimiter(temp_conserve.Energy,j,slope_limiter_type));

        //convert back to primitave variables
        solution_cell_state->Density[j] = solution_conserve.Mass[j] / dx;
        solution_cell_state->Velocity[j] = (solution_conserve.Momentum[j]/solution_cell_state->Density[j]) / dx;
        solution_cell_state->Pressure[j] = ((solution_conserve.Energy[j]/dx) - (0.5 * solution_cell_state->Density[j] * pow(solution_cell_state->Velocity[j],2.0))) * (2.0/3.0);
        
        //Conservation of Momentum = rho*u
        //Conservation of Energy = (rho * u^2) + p
        //Equation of state = u*(E + p)
    }
}

//picks a slope limiter from a list in SlopeLimiters.c
//TODO impliment slopelimiter picker

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

// C function to find maximum in arr[] of size n
double largest(cell_state temp_cell_state)
{
    int i;
    double c;
    // Initialize maximum element
    i = 0;

    c = abs(sqrt(gamma_val*temp_cell_state.Pressure[i]/temp_cell_state.Density[i]));
    double max = (abs(temp_cell_state.Velocity[i])) + c;
 
    // Traverse array elements from second and
    // compare every element with current max 
    for (i = 1; i < gridSize; i++)
        if (abs((temp_cell_state.Velocity[i]) + c) > max)
            max = abs((temp_cell_state.Velocity[i]) + c);
 
    return max;
}

//gets the variable time step dt
double getDT(cell_state temp_cell_state, double courant, double dx)
{
    double maxSignalSpeed;
    double dt;

    maxSignalSpeed = largest(temp_cell_state);
    dt = courant * (dx/maxSignalSpeed);

    printf("dt was found to be %f\n", dt);
    //system("pause");

    return dt;
}