#include <stdio.h>

int main(){
   
    double pos[] = { 0.2, 0.2, 1};
    double vel[] = { 0.0, 0.0, -1e-1};
    double cx = 0.0;
    double cy = 0.0;
    double h = 0.7;
    double r = 0.3;
    double tcoll;
    collides_with_cup(pos, vel, h, r, cx, cy, &tcoll);
    return 0;
}
