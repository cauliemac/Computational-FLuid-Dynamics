#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/*
int main ()
{
    float cents = 0.00; //current cent count
    int cycles = 10 * 1000; //in thousands
    
    for (int i = 0; i < cycles; i++)
    {
        float cents = cents + 0.01;
        printf("Was it %.2f euro?\n", cents);
        usleep(400000);
    }
}
*/

int main()
{
    sleep(3);
    for (int i = 0; i < 5; i++)
    {
        int count = 5-i;
        printf("starting gladiator games in %i\n",count);
        sleep(1);
    }
    //sleep(1);
    printf("\nAND SUPRISE SUPRISE\n");
    sleep(2);
    printf("\nTHOMAS WON");
}