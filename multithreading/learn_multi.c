#include <stdio.h>
#include <stdlib.h>
// include the open MP header file
#include <omp.h>

/* 
A program to learn the basic operations in Open MP

Thomas McCarthy
*/

//the main function
int main(int *argc, char *argv[])
{
    //specify the block to be computed in parallel
    
    #pragma omp parallel num_threads(8)
    {
        //print hello world on each thread
        printf("Hello world\n");

    }
    system("pause");
    return 0;
}