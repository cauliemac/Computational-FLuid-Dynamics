#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/*
#define gridSize 50 //size of grid
double gridSpacing = 2.0 / (gridSize);   //grid spacing ( also h)
double initialConditions[gridSize];

static int getInitialConditions(double *initialConditions, int grid, int a, int b, int sine); //gives an initial conditions array

int getInitialConditions(double *initialConditions, int gridSize, int a, int b, int sine)
{
    if (sine == 1)
        for (int i=0; i<gridSize; i++)
        {
            initialConditions[i] = sin(2*i/(gridSize/M_PI));
        }

    else
       for (int i=0; i<gridSize; i++)
        {                
            if (i/a > 1 && i/b < 1)
                initialConditions[i] = 2;
            else
                initialConditions[i] = 1;
        }
    return initialConditions
}

int main()
{
    double a = getInitialConditions(double *initialConditions, initialConditions, gridSize, 30, 60, 1);
    printf("%d",a);
    return 0;
}

/*
populates the grid between a and b (as a percent) with the height 2.
Or returns a sine wave if sine = 1
*/


#define gridsize 50
int plus = 5;
double initialConditions[gridsize];

static double add(double *initialConditions, int b);

double add(double *initialConditions ,int b)
{
    for(int i; i > gridsize; i++)
        initialConditions[i] = i+b;
    
    return *initialConditions;
}

int main()
{
    double a = add(initialConditions, plus);
    printf("%i", a);
}