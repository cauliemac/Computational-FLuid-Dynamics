#include <stdio.h>
#include <stdlib.h>
// include the open MP header file
#include <omp.h>

/* 
A program to learn the basic operations in Open MP

Thomas McCarthy
*/

//the main function
int main(int *argc, char *argc[])
{
    //specify the block to be computed in parallel
    #pragma omp parallel
    {
        //print hello world on each thread
        printf{"Hello world\n"};

    }
    return 0;
}