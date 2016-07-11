#include <stdlib.h>
#include <stdio.h>
#include <time.h>  // time()
#include <math.h>  //fabs()

#include "header.h"

int main(int argc, char **argv)
{
    struct Container container;
    container.length = LENGTH;
    container.width = WIDTH;
    container.height = HEIGHT;
    container.head = NULL;

    time_t t;
    srand((unsigned) time(&t));
    
    simulate(RB, NA, NB, &container);

    printf(test_overlap(container)? "overlap!!!\n" : "no overlap\n");
    logistic(container);

    return 0;
}

void simulate(double Rb, unsigned int Na, unsigned int Nb, struct Container *container)
{
    double r;
    unsigned int num_balls;
    double x_init, y_init;
 
    /* dropping balls 
     * Simulation terminates when either NA + NB balls are dropped 
     * or the container is full
     */
    for (num_balls = 0 ; num_balls < 10 ; num_balls++){
        /* randomly choose a ball according to NA and NB */
        r = rand_ball(NA, NB, RA, RB);
        //printf("r = %.3lf\n", r);
        /* randomly choose the initial position of the ball */
        rand_init_position(&x_init, &y_init, *container, r);
        //printf("x_init = %.3lf, y_init = %.3lf\n", x_init, y_init);
        /* construct the sphere */
        struct Sphere *sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
        sphere->radius = r;
        sphere->x = x_init;
        sphere->y = y_init;
        sphere->z = container->height;
        /* Compute the final position of the ball */
        update_sphere(sphere, container);
        /* update the container */
        update_container(sphere, container);
    }
}

void update_sphere(struct Sphere *sphere, struct Container *container)
{
    struct Sphere * first_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    struct Sphere * second_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    struct Sphere * third_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    
    first_sphere = find_first_collision(sphere, container);
    // If no other balls collide with this one, this ball will drop to the bottom
    if(first_sphere == NULL){
        sphere->z = sphere->radius;
    }
    else{
        // collision between sphere and first_sphere
        sphere->z = compute_z_collide(first_sphere, sphere);
    }
        
    
}

/*
 *  The sphere is in its initial position. Then it drops along the z axis.
 *  The function will stop the sphere the first time it collides with either
 *  another sphere or the bottom of the container.
 *  The function returns truethe pointer to the sphere that this sphere hits
 *  or NULL if the sphere hits the bottom.
 */
struct Sphere * find_first_collision(struct Sphere *sphere, struct Container *container)
{
    struct sphere_list *slist = container->head;
    struct Sphere *sphere_hit = (struct Sphere*)malloc(sizeof(struct Sphere));
    sphere_hit = NULL;
    double z_max = sphere->radius;
    double z_collide;
    while(slist != NULL){
        if (collide_z(*sphere, *(slist->sphere))){
            z_collide = compute_z_collide( slist->sphere, sphere );
            if(z_collide > z_max){
                z_max = z_collide;
                sphere_hit = slist->sphere;
            }
        }
        slist = slist->next;
    }
    return sphere_hit;
}

/*
 *  After we find the first collision, we are now trying to find the second.
 *  The first collision ball is also passed as the arg.
 *  The function returns the pointer to the found second collision sphere
 *  or NULL if no second sphere found.
 */
struct Sphere * find_second_collision(struct Sphere *sphere, struct Sphere *first_sphere, struct Container *container)
{
    
}

