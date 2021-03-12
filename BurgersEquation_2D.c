#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*
Burgers Equation using finite volume method
*/

int main ()
{   
    //create and open a file in write mode to store the initial conditian
    FILE *initial_file = NULL;
    initial_file = fopen("Burgers_2D_results/Burgers2DInitial.txt", "w");

    // declare starting variables
    //TODO set timestep size as variable!
    int gridSize = 41;   //size of grid
    double gridSpacing = 2.0 / (gridSize - 1);   //grid spacing ( aslo h)
    double evolutions = 100;  //number of evolutions
    double timestepSize = 0.01;  //size of each timestep ( also k)
    int c = 1;  //speed is 1 for simplicity

    //Create 3 arrays to store the solutions and the intitial conditions
    double* arraySolution = (double*)calloc(gridSize, sizeof(double));
    double* arrayTemp = (double*)calloc(gridSize, sizeof(double));
    double* initialConditions = (double*)calloc(gridSize, sizeof(double));

    
    /*
    Populates the first solution, i.e. the initial conditions with a height of 2 between 0.5 and 1
    everywhere else gets a height of 1
    */
   
    /* 2D array declaration*/
   double myArraySolutions_u[gridSize][gridSize];

    //TODO set centre of grid as variable

    /*set the centre of the array at 2*/
    for (int i = 0; i < gridSize - 1; i++ ) 
    {
      for (int j = 0; j < gridSize - 1; j++ ) 
      { 
        if (j >= 6 && j <= 14 && i >= 6 && i <= 14 )
          myArraySolutions_u[i][j] = 2;
        else
          myArraySolutions_u[i][j] = 1;
      }
    }
    
    for (int i = 0; i < gridSize - 1; i++ ) 
    {
      for (int j = 0; j < gridSize - 1; j++ )
      { 
        fprintf(initial_file,"myArraySolutions_u[%d][%d] = %f\n", i, j, myArraySolutions_u[i][j]);
      }
    }

    /*
    for (int i = 0; i < gridSize; i++)
    {
        //sets the initial conditions to half a wavelength
        //arraySolution[i] = sin(2 * i / (gridSize / M_PI));
        
        for (int j = 0; j < gridSize; j++)
            
            if (i * gridSpacing >= 0.4 && i * gridSpacing <= 0.6 && j * gridSpacing >= 0.4 && j * gridSpacing <= 0.6)
            {
                arraySolution[i] = 2;
            }
 
            else
            {
                arraySolution[i] = 1;
            }
    
    }
    */

    /*
    since the above for loop populates the solution array with the initial
    conditions I'm calling memcpy to copy the sol arrays contents to the
    initial condition array then writing these to the initial condition data file
    */
    double myArrayTemp_u[gridSize][gridSize];
    double myArrayTemp_v[gridSize][gridSize];
    double myArraySolutions_v[gridSize][gridSize];
    
    memcpy(initialConditions, arraySolution, gridSize * sizeof(double));

    //memcpy(myArrayTemp_u, myArraySolutions_u, sizeof(myArrayTemp_u));
    //memcpy(myArrayTemp_v, myArraySolutions_u, sizeof(myArrayTemp_u));
    memcpy(myArraySolutions_v, myArraySolutions_u, sizeof(myArrayTemp_u));
    
    /*
    for (int i = 0; i < gridSize; i++)
    {
        //print the initial conditians to file
        fprintf(initial_file, "%f\n", initialConditions[i] );
    }
    */

    //do the math that does the thing
    //loop for all evolutions
    //loop for all grids
    //calculate next value at that gridpoint
    for (int i = 0; i < evolutions; i++)
    {
        //copies the solutions array onto the temp array
        memcpy(arrayTemp, arraySolution, gridSize*sizeof(double));
        
        memcpy(myArrayTemp_u, myArraySolutions_u, sizeof(myArrayTemp_u));
        memcpy(myArrayTemp_v, myArraySolutions_v, sizeof(myArrayTemp_v));

        //changes file name with evolution cycle.
        FILE *fpointer = NULL;
        char buffer[64]; // The filename buffer.
        // Put "file" then i then ".txt" in to filename.
        snprintf(buffer, sizeof(char) * 64, "Burgers_2D_results/Burgers2DSolution%i.txt", i);
        fpointer = fopen(buffer, "w");

        //TODO fix the cfl condition here
        for (int j=1; j < gridSize-1; j++)
        {
            for (int k=1; k < gridSize-1; k++)
            {
                //Burgers 2D
                //calculate next values if u speed is pos or neg
                if (myArraySolutions_u[j][k] > 0)

                    myArraySolutions_u[j][k] = myArrayTemp_u[j][k] - timestepSize * 
                    
                    ((myArrayTemp_u[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_u[j][k],2)) - (0.5 * pow(myArrayTemp_u[j-1][k],2))) +

                    ((myArrayTemp_v[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_u[j][k],2)) - (0.5 * pow(myArrayTemp_u[j][k-1],2)))));

                else if (myArraySolutions_u[j][k] < 0)

                    myArraySolutions_u[j][k] = myArrayTemp_u[j][k] - timestepSize * 
                    
                    ((myArrayTemp_u[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_u[j+1][k],2)) - (0.5 * pow(myArrayTemp_u[j][k],2))) +

                    ((myArrayTemp_v[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_u[j][k+1],2)) - (0.5 * pow(myArrayTemp_u[j][k],2)))));
                
                //calculate next values if v speed is pos or neg
                else if (myArraySolutions_v > 0)

                    myArraySolutions_v[j][k] = myArrayTemp_v[j][k] - timestepSize * 
                    
                    ((myArrayTemp_u[j][k] / gridSpacing)*(((0.5 * pow(myArrayTemp_v[j][k],2)) - (0.5 * pow(myArrayTemp_v[j-1][k],2))) +

                    ((myArrayTemp_v[j][k] / gridSpacing)*((0.5* pow(myArrayTemp_v[j][k],2)) - (0.5 * pow(myArrayTemp_v[j][k-1],2))))));

                else if (myArraySolutions_v < 0)

                    myArraySolutions_v[j][k] = myArrayTemp_v[j][k] - timestepSize * 
                    
                    ((myArrayTemp_u[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_v[j+1][k],2)) - (0.5 * pow(myArrayTemp_v[j][k],2))) +

                    ((myArrayTemp_v[j][k] / gridSpacing)*((0.5 * pow(myArrayTemp_v[j][k+1],2)) - (0.5 * pow(myArrayTemp_v[j][k],2)))));
                

                //else
                    //myArraySolutions_u[j][k] = myArrayTemp_u[j][k] - (timestepSize / gridSpacing) * 
                    
                    //(((0.5 * pow(myArrayTemp_u[j+1][k+1], 2)) - (0.5 * pow(myArrayTemp_u[j][k], 2))));


                 //print the x axis label (which is j) and the solution to a text file
                 fprintf(fpointer,"%i \t %i \t %f \n", j, k, myArraySolutions_u[j][k]);
            }
        }      
    }
    printf("grid Spacing: %f\n", gridSpacing);
    printf("evolutions: %f\n", evolutions);
    printf("time step size: %f\n", timestepSize);

    return EXIT_SUCCESS;
}