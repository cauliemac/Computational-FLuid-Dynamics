#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SlopeLimiters.h"

/*
1D Linear Advection using finite volume method
With implimentation of slope limiters


TODO move each section to it's own function
TODO add slopelimiters to header files

*/

double slopeLimiter(double sigma);    //declare the slope limiter function
double slopeLimiter_L(double sigma);


int main ()
{   
    //create and open a file in write mode to store the initial conditian
    FILE *initial_file = NULL;
    initial_file = fopen("LinearAdvection_1D_results/LinearAdvectionInitial.txt", "w");

    // declare starting variables
    int gridSize = 50;   //size of grid
    double gridSpacing = 2.0 / (gridSize );   //grid spacing ( also h)
    double evolutions = 100;  //number of evolutions
    double timestepSize = 0.03;  //size of each timestep ( also k)
    int c = 1; //speed of equation 1 for simple
    double courant = timestepSize/ gridSpacing; //Cournant number for printout
    double minmod_a;
    double minmod_b;
    double minmod_b_L;
    double minmod_a_L;
    double sigma;
    double sigma_L;

    //Create 3 arrays to store the solutions and the intitial conditions
    double arraySolution[gridSize];
    double arrayTemp[gridSize];
    double initialConditions[gridSize];
    
    /*
    Populates the grid with initial conditions of height 2
    everywhere else gets a height of 1
    or populates the grid with a sin curve
    */
    for (int i = 0; i < gridSize; i++)
    {
        //sets the initial conditions to one wavelength of sin
        //arraySolution[i] = sin(2 * i / (gridSize / M_PI));
                
        if (i*gridSpacing >= 0.3 && i*gridSpacing <= 0.9)
            arraySolution[i] = 2;
        else
            arraySolution[i] = 1;
    }

    //memcpy to copy initial conditioan onto the solution array
    memcpy(initialConditions, arraySolution, gridSize * sizeof(double));

    for (int i = 0; i < gridSize; i++)
    {
        //print the initial conditians to file
        fprintf(initial_file, "%f\n", initialConditions[i] );
    }


    //do the math that does the thing
    //loop for all evolutions
    //loop for all grid points
    //calculate next value at that gridpoint
    for (int i = 0; i < evolutions; i++)
    {
        //copies the solutions array onto the temp array
        memcpy(arrayTemp, arraySolution, gridSize*sizeof(double));

        //changes file name with evolution cycle.
        FILE *fpointer = NULL;
        char buffer[64]; // The filename buffer.

        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 64, "LinearAdvection_1D_results/LinearAdvectionSolution%i.txt", i);
        fpointer = fopen(buffer, "w");

        //calculates the next value of the current cell
        for (int j = 2; j < gridSize-1; j++)
        {
            //For the size of the grid, calculate the 2nd order of the Linear Advection Equation
            if (arraySolution[j] > 0)
                //uses the minmod slope limiter to calculate the values for sigma and sigma-1
                minmod_a = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing;
                minmod_b = (arrayTemp[j+1] - arrayTemp[j])/gridSpacing;

                if(minmod_a * minmod_b < 0)
                    sigma = 0;
                else if (abs(minmod_a) < abs(minmod_b) && minmod_a * minmod_b > 0)
                    sigma = minmod_a;
                else if (abs(minmod_b) < abs(minmod_a) && minmod_a * minmod_b > 0)
                    sigma = minmod_b;

                minmod_a_L = (arrayTemp[j-1] - arrayTemp[j-2])/gridSpacing;
                minmod_b_L = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing;

                if(minmod_a_L * minmod_b_L < 0)
                    sigma_L = 0;
                else if (abs(minmod_a_L) < abs(minmod_b_L) && minmod_a_L * minmod_b_L > 0)
                    sigma_L = minmod_a_L;
                else if (abs(minmod_b_L) < abs(minmod_a_L) && minmod_a_L * minmod_b_L > 0)
                    sigma_L = minmod_b_L;

                // calculates next value
                //arraySolution[j] = arrayTemp[j] - (timestepSize/gridSpacing) * ((0.5 * pow(arrayTemp[j], 2) - (0.5 * pow(arrayTemp[j-1], 2))));
                arraySolution[j] = arrayTemp[j] - courant * c * (arrayTemp[j] - arrayTemp[j-1]) - 0.5 * c * courant * (gridSpacing - timestepSize)*(sigma - sigma_L);

            //for negative velocity
            //TODO this needs fixed. Else function broken
            /*
            else
                //next value
                arraySolution[j] = arrayTemp[j] - (timestepSize / gridSpacing) * (((0.5 * pow(arrayTemp[j + 1], 2)) - (0.5 * pow(arrayTemp[j], 2))));
            */
            
            //uncomment below to print values to console
            //printf("%f\n", arraySolution[j])
            
            //print the x axis label (which is j) and the solution to a text file
            fprintf(fpointer, "%i \t %f\n", j, arraySolution[j]);
        }
    }
    printf("grid Spacing: %f\n", gridSpacing);
    printf("time step size: %f\n", timestepSize);
    printf("evolutions: %f\n", evolutions);
    printf("Courant number: %f\n", courant);

    return EXIT_SUCCESS;
}

/*
//GET THESE WORKING TO TIDY UP CODE
double slopeLimiter(double sigma)   //minmod
    {
        double a;
        double b;

        a = (arrayTemp[j] - arrayTemp[j-1])/gridSpacing);
        b = (arrayTemp[j+1] - arrayTemp[j])/gridSpacing);
        
        result_Q = ;
        
        
        if (result_Q > fmax(a,b))
            result_Q = fmax(a,b);
        
        return result_Q;
    }

double slopeLimiter_L(double sigma)
    {
        double result_Q;
        
        double theta = ratio(a,b);



        result_Q = ;
        
        
        if (result_Q > fmax(a,b))
            result_Q = fmax(a,b);
        
        return result_Q;
    }
*/

//returns the ratio of the slopes from the cells
//to the right and left of a cell
/*
double ratio(double cell_L, double cell_R)
    {
        double ratio;

        ratio = (cell_R - cell_L)/(cell_L - cell_R);

        return ratio;
    }
*/