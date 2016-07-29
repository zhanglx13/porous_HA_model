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
    
    logistic(container);

    printf(test_overlap(container)? "overlap!!!\n" : "no overlap\n");

    /* struct Sphere s0, s1; */
    /* s0.x = 50.0; */
    /* s0.y = 80.0; */
    /* s0.z = 110.0; */
    /* s0.radius = 45; */
    /* s0.id = 0; */
    /* s1.x = 150; */
    /* s1.y = 90.0; */
    /* s1.z = 40; */
    /* s1.radius = 35; */
    /* s1.id = 1; */
    
    
    /* double x, y, r = 27.0; */
    /* double z = z_collide_two(&s0, &s1, r, &x, &y); */
    /* if (z == -1) */
    /*     return 0; */
    /* double err1 = (x-s0.x)*(x-s0.x)+(y-s0.y)*(y-s0.y)+(z-s0.z)*(z-s0.z) - (r + s0.radius)*(r + s0.radius); */
    /* double err2 = (x-s1.x)*(x-s1.x)+(y-s1.y)*(y-s1.y)+(z-s1.z)*(z-s1.z) - (r + s1.radius)*(r + s1.radius); */
    
    /* printf("err1 = %.3lf\nerr2 = %.3lf\n", err1, err2); */

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
    for (num_balls = 0 ; num_balls < NBALLS ; num_balls++){

        //printf("Ball %u is coming ... \n", num_balls);
        /* randomly choose a ball according to NA and NB */
        r = rand_ball(NA, NB, RA, RB);
        //r = RB;
        //printf("r = %.3lf\n", r);
        /* randomly choose the initial position of the ball */
        rand_init_position(&x_init, &y_init, *container, r);
        //x_init = 500.0;
        //y_init = 500.0;
        //printf("x_init = %.3lf, y_init = %.3lf\n", x_init, y_init);
        /* construct the sphere */
        struct Sphere *sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
        sphere->radius = r;
        sphere->x = x_init;
        sphere->y = y_init;
        sphere->z = container->height;
        sphere->id = num_balls;
        /* Compute the final position of the ball */
        update_sphere(sphere, container);
        /* update the container */
        update_container(sphere, container);
    }

    /* num_balls = 1; */
    /*     x_init = 510.0; */
    /*     y_init = 500.0; */
    /*     r = RA; */
    /*     struct Sphere *sphere1 = (struct Sphere*)malloc(sizeof(struct Sphere)); */
    /*     sphere1->radius = r; */
    /*     sphere1->x = x_init; */
    /*     sphere1->y = y_init; */
    /*     sphere1->z = container->height; */
    /*     sphere1->id = num_balls; */
    /*     /\* Compute the final position of the ball *\/ */
    /*     update_sphere(sphere1, container); */
    /*     /\* update the container *\/ */
    /*     update_container(sphere1, container); */
        
        


      
}

void update_sphere(struct Sphere *sphere, struct Container *container)
{
    struct Sphere * first_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    
    //struct Sphere * third_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    
    first_sphere = find_first_collision(sphere, container);
    // If no other balls collide with this one, this ball will drop to the bottom
    if(first_sphere == NULL){
        sphere->z = sphere->radius;
        free(first_sphere);
    }
    else{
        // collision between sphere and first_sphere
        sphere->z = compute_z_collide(first_sphere, sphere);
        if (overlap(*first_sphere, *sphere))
            printf("First_sphere %u overlaps with sphere %u\n", first_sphere->id, sphere->id);
        printf("ball %u(%.3lf, %.3lf, %.3lf, %.3lf) and %u(%.3lf, %.3lf, %.3lf, %.3lf) collided\n", 
               first_sphere->id, first_sphere->x, first_sphere->y, first_sphere->z, first_sphere->radius, 
               sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
        struct Sphere * second_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
        int flag = -1;
        second_sphere = find_second_collision(sphere, first_sphere, container, &flag);

        if (second_sphere == NULL){
            if (flag == -1){
                printf("flag == -1\n");
                update_sphere(sphere, container);
            }
            else if(flag == HIT_BOTTOM){
                printf("Hit bottom\n");
                return ;
            }
            else{
                printf("flag == hit_wall\n");
                find_third_collision();
            }
                
        }else{
            printf("hit second_sphere\n");
            find_third_collision();
        }

    }
        
    
}

struct Sphere *find_third_collision()
{
    return NULL;
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
struct Sphere * find_second_collision(struct Sphere *sphere, struct Sphere *first_sphere, struct Container *container, int *flag)
{
    /****************************************************************
     * stage 1:
     *     Find a free position of the sphere around first_sphere
     *     by traversing all the positions with delta_x.
     *     If found, return NULL. If cannot find one, go to stage 2.
     ***************************************************************/
    double r0 = first_sphere->radius;
    double r1 = sphere->radius;
    double start_pos = first_sphere->x + r0 + r1;
    double end_pos = first_sphere->x - r0 - r1;
    double y_cord;
    // save the initial position of sphere
    double x_init = sphere->x;
    double y_init = sphere->y;
    double z_init = sphere->z;
    // initial position
    sphere->x = start_pos; // x1 = x0 + r0 + r1
    sphere->y = first_sphere->y;
    sphere->z = first_sphere->z;
    // step setup
    double steps = 1000.0;  // number of steps to scan half cycle
    double range = 2.0 * (r0 + r1);
    double delta_x = range / steps;
    // test initial position
    if (test_sphere_safety(container, sphere))
            return NULL;
    sphere->x = sphere->x - delta_x;
    while(sphere->x >= end_pos){
        y_cord = sqrt( (r0 + r1 )*(r0 + r1) - (sphere->x - first_sphere->x)*(sphere->x - first_sphere->x) );
        // test the up position
        sphere->y = first_sphere->y + y_cord;
        if (test_sphere_safety(container, sphere))
            return NULL;
        // test the down position
        sphere->y = first_sphere->y - y_cord;
        if (test_sphere_safety(container, sphere))
            return NULL;
        // if the current position is not safe
        // compute the next position
        sphere->x = sphere->x - delta_x;
    }

    /********************************************************************
     * stage 2:
     *     If no free position is found in stage 1, we try to test if 
     *     the shere can hit the bottom
     ********************************************************************/
    // When the sphere can hit the bottom, it must satisfy the following conditions:
    // 1. sphere->radius  >  first_sphere->radius
    // 2. first_sphere->z = first_sphere->radius
    if ( (r1 > r0) && (first_sphere->z == r0) ){
        double r = sqrt( (r0+r1)*(r0+r1) - (r1-r0)*(r1-r0) );
        start_pos = first_sphere->x + r;
        end_pos = first_sphere->x -r;
        // initial position
        sphere->x = start_pos;
        sphere->y = first_sphere->y;
        sphere->z = r1;
        // step setup
        steps = 1000.0;
        range = 2.0 * r;
        delta_x = range / steps;
        // test initial position
        if (test_sphere_safety(container, sphere)){
            *flag = HIT_BOTTOM;
            return NULL;
        }
        sphere->x = sphere->x - delta_x;
        while(sphere->x >= end_pos){
            y_cord = sqrt( r*r - (sphere->x - first_sphere->x)*(sphere->x - first_sphere->x) );
            // test the up position
            sphere->y = first_sphere->y + y_cord;
            if (test_sphere_safety(container, sphere)){
                *flag = HIT_BOTTOM;
                return NULL;
            }
            // test the down position
            sphere->y = first_sphere->y - y_cord;
            if (test_sphere_safety(container, sphere)){
                *flag = HIT_BOTTOM;
                return NULL;
            }
            // if the current position is not safe
            // compute the next position
            sphere->x = sphere->x - delta_x;
        }
    }

    /********************************************************************
     * stage 3:
     *     If the sphere cannot hit the bottom without colliding any
     *     other spheres, we need to find the second collision
     ********************************************************************/
    printf("no available space\n");
    // restore the initial position of the sphere
    sphere->x = x_init;
    sphere->y = y_init;
    sphere->z = z_init;
    
    struct sphere_list *slist = container->head;
    struct Sphere *sphere_hit = (struct Sphere*)malloc(sizeof(struct Sphere));
    sphere_hit = NULL;
    double z_min = sphere->z;
    double z_collide;
    double dist;
    double x = sphere->x, y = sphere->y;
    double x_temp, y_temp;
    while( slist != NULL ){
        if (slist->sphere->id != first_sphere->id){
            dist = (slist->sphere->x - first_sphere->x)*(slist->sphere->x - first_sphere->x) +
                   (slist->sphere->y - first_sphere->y)*(slist->sphere->y - first_sphere->y) +
                   (slist->sphere->z - first_sphere->z)*(slist->sphere->z - first_sphere->z);
            // only test those spheres that might collide with the sphere
            if ( dist <= (r0 + r1 + r1)*(r0 + r1 + r1) ){
                z_collide = z_collide_two(first_sphere, slist->sphere, r1, &x_temp, &y_temp);
                if (z_collide < z_min){
                    z_min = z_collide;
                    x = x_temp;
                    y = y_temp;
                    sphere_hit = slist->sphere;
                }
            }
        }
        slist = slist->next;
    }
    // The new position of the sphere will be the lowest one along z axis
    sphere->x = x;
    sphere->y = y;
    sphere->z = z_min;
    
    //
    // Here we consider if the sphere can hit the wall
    //
    if ( (first_sphere->x + r0 + r1 + r1 > container->length)  || (first_sphere->y + r0 + r1 + r1 > container->width) ){
        // Hit the x limit
        if (first_sphere->x + r0 + r1 + r1 > container->length){
            double AE = container->length - first_sphere->x - r1;
            double BE = sqrt( (r0+r1)*(r0+r1) - AE*AE );
            z_collide = first_sphere->z + BE;
            if (z_collide < z_min){
                z_min = z_collide;
                x = first_sphere->x + AE;
                y = first_sphere->y;
                *flag == HIT_WALL_X;
            }
        }

        // Hit the y limit
        if (first_sphere->y + r0 + r1 + r1 > container->width){
            double AE = container->width - first_sphere->y - r1;
            double BE = sqrt( (r0+r1)*(r0+r1) - AE*AE );
            z_collide = first_sphere->z + BE;
            if (z_collide < z_min){
                z_min = z_collide;
                x = first_sphere->x;
                y = first_sphere->y + AE;
                *flag == HIT_WALL_Y;
            }
        }
        
        if (*flag == HIT_WALL_X  || *flag == HIT_WALL_Y){
            sphere->x = x;
            sphere->y = y;
            sphere->z = z_min;
            return NULL;
        }
    }
    
    //
    // sphere_hit should be asserted != NULL
    //
    if (sphere_hit == NULL)
        printf("*************sphere_hit = NULL!!!!**************\n");
    
    return sphere_hit;
}




