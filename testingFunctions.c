#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Windows.h>

#define gridSize 50

typedef struct CELL_STATE
{
    double Density[gridSize];
    double Pressure[gridSize];
    double Velocity[gridSize];
}cell_state;

void set(cell_state a);//declare

void add(cell_state a, cell_state b, cell_state *result);//declare

int main()
{
    cell_state temp, solution, result;

    for(int i=1; i<gridSize+1;i++)
    {
        temp.Density[i] = i;
        temp.Pressure[i] = i;
        temp.Velocity[i] = i;

        printf("set temp.Density[%i] = %d\n", i, temp.Density[i]);
        Sleep(20);
    }

    //set(temp);
    //set(solution);

    add(temp, solution, &result);

    system("pause");

    return 0;
}

/*
void set(cell_state a)
{
    for(int i=0; i<=gridSize;i++)
    {
        a.Density[i] = i;
        a.Pressure[i] = i;
        a.Velocity[i] = i;

        printf("set a.Density[%i] = %d\n", i, &a.Density[i]);
        Sleep(20);
    }
    //printf("%d\n", a.Density[gridSize]);
    Sleep(2000);
}
*/

void add(cell_state a, cell_state b, cell_state *result)
{
    result->Density[gridSize] = a.Density[gridSize] + b.Density[gridSize];
    printf("result->Density[%i] = %d\n", gridSize, result->Density[gridSize]);
    Sleep(2000);
}