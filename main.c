#include <stdlib.h>
#include <stdio.h>
#include <time.h>  // time()
#include <math.h>  //fabs()

#include "header.h"

#define RAND0(bound)        \
    (double)rand()*(double)(bound) / (double)RAND_MAX
    


#define LENGTH   150
#define HEIGHT   150
#define WIDTH    150

#define RA       5
#define RB       50

#define NA      50
#define NB      100



int main(int argc, char **argv)
{
    struct Container container;
    container.length = LENGTH;
    container.width = WIDTH;
    container.height = HEIGHT;
    container.head = NULL;

    time_t t;
    srand((unsigned) time(&t));
    
    

    simulate(RB, NA, NB, container);


    struct Sphere s0, s1;
    s0.x = 5;
    s0.y = 5;
    s0.z = 0;
    s0.radius = 2;
    s1.radius = 3;
    s1.x = 8;
    s1.y = 9;
    s1.z = 0;
    printf(collide_z(s0, s1) ? "collide\n" : "fine\n");


    printf("%d\n", LENGTH*WIDTH*HEIGHT);

    return 0;
}

void simulate(double Rb, unsigned int Na, unsigned int Nb, struct Container container)
{
    double r;
    unsigned int num_balls;
    double x_init, y_init;
    double x_final, y_final;
    /* dropping balls 
     * Simulation terminates when either NA + NB balls are dropped 
     * or the container is full
     */
    for (num_balls = 0 ; num_balls < NA+NB ; num_balls++){
        /* randomly choose a ball according to NA and NB */
        r = rand_ball(NA, NB, RA, RB);
        /* randomly choose the initial position of the ball */
        rand_init_position(&x_init, &y_init, container, r);
        /* construct the sphere */
        struct Sphere sphere;
        sphere.radius = r;
        sphere.x = x_init;
        sphere.y = y_init;
        sphere.z = container.height;
        /* Compute the final position of the ball */
        update_sphere(&sphere, container);
        /* update the container */
        update_container(sphere, &container);
    }
}


/*
  Pick a ball randomly. The probability is proportional to 
  the number of balls. The function returns the radius of 
  the picked ball
  @param Na: number of ball A
  @param Nb: number of ball B
  @param Ra: radius of ball A
  @param Rb: radius of ball B
  @return radius of the picked ball
 */
double rand_ball(unsigned int Na, unsigned int Nb, double Ra, double Rb)
{
    double rd = RAND0(Na + Nb);
    if (rd <= Na)
        return Ra;
    else
        return Rb;
}

/*
 *  Choose a random initial position for the ball
 */
void rand_init_position(double *x, double *y, struct Container container, double r)
{
    *x = RAND0(container.length - 2*r);
    *x += r;
    *y = RAND0(container.width - 2*r);
    *y += r;
}

void update_sphere(struct Sphere *sphere, struct Container container)
{
    
}

void update_container(struct Sphere sphere, struct Container *container)
{
    
}

/*
 *   Test whether two spheres collide when projected in the x-y plane
 */
bool collide_z(struct Sphere s0, struct Sphere s1)
{
    double delta_x = s0.x - s1.x;
    double delta_y = s0.y - s1.y;
    return sqrt(delta_x * delta_x + delta_y * delta_y) < (s0.radius + s1.radius) ? true : false; 
}
