/* C graphics program to draw a line */
#include <stdlib.h>
#include <stdio.h>
#include <graphics.h>
#include <conio.h>

int main() {
    int gd = DETECT, gm;
    /* initialization of graphic mode */
    initgraph(&gd, &gm, "C:\\TC\\BGI"); 
    line(100,100,200, 200);
    getch();
    closegraph();
    return 0;
}