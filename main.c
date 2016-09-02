#include <stdlib.h>
#include <stdio.h>
#include <time.h>  // time()
#include <math.h>  //fabs()

#include "header.h"

int main(int argc, char **argv)
{
    

    time_t t;
    srand((unsigned) time(&t));
    
    run(5.1);
    /* run(10); */
    /* run(20); */
    /* run(30); */
    /* run(50); */
    
    
    return 0;
}

void simulate(double Rb, unsigned int Na, unsigned int Nb)
{
    double r;
    unsigned int num_balls, ind = 0, num_A = 0, num_B = 0;
    double x_init, y_init;
    struct Container *container = (struct Container*)malloc(sizeof(struct Container));
    container->length = LENGTH;
    container->width = WIDTH;
    container->height = HEIGHT;

    //    do {
    
        container->head = NULL;
        num_A = 0;
        num_B = 0;
        ind = 0;
        x_init = 1.0;
        y_init = 1.0;
   
        /* dropping balls 
         * Simulation terminates when either NA + NB balls are dropped 
         * or the container is full
         */
        for (num_balls = 0 ; num_balls < NBALLS ; num_balls ++){
            /* randomly choose a ball according to NA and NB */
            r = rand_ball(Na, Nb, RA, Rb);
            /* randomly choose the initial position of the ball */
            rand_init_position(&x_init, &y_init, *container, r);

            /* x_init += 2.0*r; */
            /* if (x_init > container->length - r){ */
            /*     x_init = r + 1.0; */
            /*     y_init += r*2.0; */
            /*     if (y_init > container->width - r) */
            /*         y_init = r + 1.0; */
            /* } */


            /* construct the sphere */
            struct Sphere *sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
            sphere->radius = r;
            sphere->x = x_init;
            sphere->y = y_init;
            sphere->z = container->height;
            sphere->id = ind;
            
            printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) is comming\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
            
            /* Compute the final position of the ball */
            update_sphere(sphere, container);
            /* Do not add the ball if the container is full */
            if (sphere->z + sphere->radius > container->height){
                //printf("Ball %u not added!\n", num_balls);
                continue;
            }
            /* update the container */
            update_container(sphere, container);
            ind ++;
            if (r == RA  && ( sphere->z <= REAL_H) )
                num_A ++;
            else if (r == Rb  &&  (sphere->z <= REAL_H))
                num_B ++;
            
            /* if (test_overlap(*container)) */
            /*     exit(1); */
            //else 
                /* printf("ball %u: (%.3lf, %.3lf, %.3lf, %.3lf)\n",  */
                /*        sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius); */
        }   

        //    }while (test_overlap(*container));
    printf("%u, %u\t", num_A, num_B);
    printf("Nb/Na = %.5lf ---> ", (double)num_B / (double)num_A);
    double Vb = (double)num_B * 4.0/3.0 * PI * Rb*Rb*Rb;
    double V = container->length * container->width * REAL_H;
    printf("%.3lf, %.0lf\t", Vb, V);
    printf("Vb/V = %.5lf\n", Vb / V);
    
}

void update_sphere(struct Sphere *sphere, struct Container *container)
{
    struct Sphere * first_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
    
    first_sphere = find_first_collision(sphere, container);
    // If no other balls collide with this one, this ball will drop to the bottom
    if(first_sphere == NULL){
        sphere->z = sphere->radius;
        free(first_sphere);
        if (!test_sphere_safety(container, sphere)){
            printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) not safe to added to the container at the end of find_first\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
        }
        else 
            printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) is safe to added to the container at the end of find_first\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
    }
    else{
        if (!test_sphere_safety(container, sphere)){
            printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) not safe to added to the container at the end of find_first\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
        }
        else 
            printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) is safe to added to the container at the end of find_first\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
        // collision between sphere and first_sphere
        sphere->z = compute_z_collide(first_sphere, sphere);
        //printf("Finding second collision for sphere %u (%.3lf, %.3lf, %.3lf, %.3lf)\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
        if (overlap(*first_sphere, *sphere))
            printf("First_sphere %u overlaps with sphere %u\n", first_sphere->id, sphere->id);
        /* printf("ball %u(%.3lf, %.3lf, %.3lf, %.3lf) and %u(%.3lf, %.3lf, %.3lf, %.3lf) collided\n",  */
        /*        first_sphere->id, first_sphere->x, first_sphere->y, first_sphere->z, first_sphere->radius,  */
        /*        sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius); */
        struct Sphere * second_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
        int flag = -1;
        second_sphere = find_second_collision(sphere, first_sphere, container, &flag);

        if (second_sphere == NULL){
            if (flag == -1){
                //printf("flag == -1\n");
                update_sphere(sphere, container);
            }
            else if(flag == HIT_BOTTOM){
                //printf("Hit bottom\n");
                return ;
            }
            else{ // HIT_WALL
                //printf("flag == hit_wall\n");
                //find_third_collision();
            }
                
        }else{
            //printf("hit second_sphere %u\n", second_sphere->id);
            //struct Sphere * third_sphere = (struct Sphere*)malloc(sizeof(struct Sphere));
            //find_third_collision(sphere, first_sphere, second_sphere, container);
        }

    }
        
    
}

bool find_second_collision_stage_1(struct Sphere *sphere, struct Sphere *first_sphere, struct Container *container)
{
    if (test_overlap(*container)){
        printf("Overlap detected at the beginning of stage 1 for sphere %u\n", sphere->id);
        exit(1);
    }
    if (!test_sphere_safety(container, sphere)){
        printf("Sphere %u not safe to add to the container at stage 1\n", sphere->id);
        exit(0);
    }
    bool foundFreeSpace = false;
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
    if (test_sphere_safety(container, sphere)){
        //printf("Found a free space!!!\n");
        return true;
    }
    sphere->x = sphere->x - delta_x;
    while(sphere->x >= end_pos){
        y_cord = sqrt( (r0 + r1 )*(r0 + r1) - (sphere->x - first_sphere->x)*(sphere->x - first_sphere->x) );
        // test the up position
        sphere->y = first_sphere->y + y_cord;
        if (test_sphere_safety(container, sphere)){
            //printf("Found a free space!!!\n");
            return true;
        }
        // test the down position
        sphere->y = first_sphere->y - y_cord;
        if (test_sphere_safety(container, sphere)){
            //printf("Found a free space!!!\n");
            return  true;
        }
        // if the current position is not safe
        // compute the next position
        sphere->x = sphere->x - delta_x;
    }
    if (foundFreeSpace == false){
        sphere->x = x_init;
        sphere->y = y_init;
        sphere->z = z_init;
    }
    if (test_overlap(*container)){
        printf("Overlap detected at the end of stage 1 for sphere %u\n", sphere->id);
        exit(1);
    }
    return foundFreeSpace;
}

bool find_second_collision_stage_2(struct Sphere *sphere, struct Sphere *first_sphere, struct Container *container, int *flag)
{
    if (test_overlap(*container)){
        printf("Overlap detected at the beginning of stage 2 for sphere %u\n", sphere->id);
        exit(1);
    }
    if (!test_sphere_safety(container, sphere)){
        printf("Sphere %u not safe to add to the container at stage 2\n", sphere->id);
        exit(0);
    }
    bool hitBottom = false;
    double r0 = first_sphere->radius;
    double r1 = sphere->radius;
    double start_pos, end_pos, y_cord;
    double steps, range, delta_x;
    // save the initial position of sphere
    double x_init = sphere->x;
    double y_init = sphere->y;
    double z_init = sphere->z;
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
            //printf("Hit the bottom!!!\n");
            return true;
        }
        sphere->x = sphere->x - delta_x;
        while(sphere->x >= end_pos){
            y_cord = sqrt( r*r - (sphere->x - first_sphere->x)*(sphere->x - first_sphere->x) );
            // test the up position
            sphere->y = first_sphere->y + y_cord;
            if (test_sphere_safety(container, sphere)){
                *flag = HIT_BOTTOM;
                //printf("Hit the bottom!!!\n");
                return true;
            }
            // test the down position
            sphere->y = first_sphere->y - y_cord;
            if (test_sphere_safety(container, sphere)){
                *flag = HIT_BOTTOM;
                //printf("Hit the bottom!!!\n");
                return true;
            }
            // if the current position is not safe
            // compute the next position
            sphere->x = sphere->x - delta_x;
        }
    }
    if (hitBottom == false){
        sphere->x = x_init;
        sphere->y = y_init;
        sphere->z = z_init;
    }
    if (test_overlap(*container)){
        printf("Overlap detected at the end of stage 2 for sphere %u\n", sphere->id);
        exit(1);
    }
    return hitBottom;
}

struct Sphere *find_second_collision_stage_3(struct Sphere *sphere, struct Sphere *first_sphere, struct Container *container, int *flag)
{
    //printf("no available space\n");
    // restore the initial position of the sphere
    /* sphere->x = x_init; */
    /* sphere->y = y_init; */
    /* sphere->z = z_init; */
    
    double steps, range;
    double r0 = first_sphere->radius;
    double r1 = sphere->radius;

    if (!test_sphere_safety(container, sphere)){
        printf("Sphere %u not safe to add to the container at stage 3\n", sphere->id);
        exit(0);
    }
    

    double x0 = first_sphere->x;
    double y0 = first_sphere->y;
    double z0 = first_sphere->z;
    
    double x1 = sphere->x;
    double y1 = sphere->y;
    double z1 = sphere->z; 

    // In the alpha frame
    double w = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
    double yA = z1-z0;
    // step setup
    steps = 1000.0;
    range = yA;
    double delta_y = range / steps;
    double xA;
    struct Sphere *sphere_hit = (struct Sphere*)malloc(sizeof(struct Sphere));
    while (yA >= 0){
        if ( ((r0+r1)*(r0+r1) - yA*yA)*((r0+r1)*(r0+r1) - yA*yA) < 0.00001  )
            xA = 0;
        else
            xA = sqrt( (r0+r1)*(r0+r1) - yA*yA );
        sphere->z = z0 + yA;
        if (w*w < 0.00001){ // w == 0
            sphere->x = x0 + xA;
            sphere->y = y0;
        }else {
            sphere->x = xA*(x1-x0) / w + x0;
            sphere->y = xA*(y1-y0) / w + y0;
        }
        sphere_hit = first_collision(sphere, container, flag);
        if (sphere_hit != NULL){
            //printf("Found second collision ball with %u\n", sphere_hit->id);
            if (test_overlap(*container)){
                printf("Overlap detected at the end of stage 3 for sphere %u\n", sphere->id);
                exit(1);
            }
            return sphere_hit;
        }
        if (*flag == HIT_WALL){
            printf("ball %u (%.3lf, %.3lf, %.3lf, %.3lf) hit the wall\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
            printf("xA = %.3lf\n", xA);
            printf("yA = %.3lf\n", yA);
            printf("r0 = %.3lf\n", r0);
            printf("r1 = %.3lf\n", r1);
            printf("w = %.3lf\n", w);
            printf("x0 = %.3lf\n", x0);
            printf("y0 = %.3lf\n", y0);
            printf("z0 = %.3lf\n", z0);
            printf("x1 = %.3lf\n", x0);
            printf("y1 = %.3lf\n", y1);
            printf("z1 = %.3lf\n", z1);
            printf("first_sphere:  %u (%.3lf, %.3lf, %.3lf, %.3lf)\n", first_sphere->id, first_sphere->x, first_sphere->y, first_sphere->z, first_sphere->radius);
            
            if (test_overlap(*container)){
                printf("Overlap detected at the end of stage 3 for sphere %u\n", sphere->id);
                exit(1);
            }

            return NULL;
        }
        yA = yA - delta_y;
    }

    if (test_overlap(*container)){
        printf("Overlap detected at the end of stage 3 for sphere %u\n", sphere->id);
        exit(1);
    }
    if (!test_sphere_safety(container, sphere)){
        printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) not safe to added to the container at the end of stage 3\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
    }
    else 
        printf("Sphere %u (%.3lf, %.3lf, %.3lf, %.3lf) is safe to added to the container at the end of stage 3\n", sphere->id, sphere->x, sphere->y, sphere->z, sphere->radius);
    
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
    //printf("Start find second collision!!!\n");
    
    /****************************************************************
     * stage 1:
     *     Find a free position of the sphere around first_sphere
     *     by traversing all the positions with delta_x.
     *     If found, return NULL. If cannot find one, go to stage 2.
     ***************************************************************/
    if (find_second_collision_stage_1(sphere, first_sphere, container))
        return NULL;
    
    /********************************************************************
     * stage 2:
     *     If no free position is found in stage 1, we try to test if 
     *     the shere can hit the bottom
     ********************************************************************/
    if (find_second_collision_stage_2(sphere, first_sphere, container, flag))
        return NULL;

    /********************************************************************
     * stage 3:
     *     If the sphere cannot hit the bottom without colliding any
     *     other spheres, we need to find the second collision
     ********************************************************************/
    return find_second_collision_stage_3(sphere, first_sphere, container, flag);
    
}
    
/* struct Sphere *find_third_collision(struct Sphere *sphere, struct Sphere *first_sphere, struct Sphere *second_sphere, struct Container *container) */
/* { */
/*     printf("Start finding the third collision!!!\n"); */

/*     double z_max = compute_z_max(sphere, first_sphere, second_sphere);; */
/*     double z_min = sphere->radius; */
/*     // step setup */
/*     double steps = 100000.0; */
/*     double range = z_max - z_min; */
/*     double delta_z = range / steps; */
/*     // setup flags */
/*     int state = 0; */
/*     // initial position */
/*     sphere->z = z_max; */
/*     double x_temp, y_temp; */
/*     double z_collide = z_collide_two(first_sphere, second_sphere, sphere->radius, &x_temp, &y_temp); */
/*     if (z_collide == -1) */
/*         printf("Not a triangle!!!\n"); */
/*     sphere->x = x_temp; */
/*     sphere->y = y_temp; */
/*     // test initial position */
/*     if (!test_sphere_safety(container, sphere)){ */
/*         printf("Initial position not safe!!!\n"); */
/*         state = 1; */
/*     } */
/*     sphere->z = sphere->z - delta_z; */
/*     double x1,x2,y1,y2; */
/*     // test left position */
/*     while (sphere->z >= z_min){ */
/*         solve_x_y(sphere, first_sphere, second_sphere, sphere->z, &x1, &y1, &x2, &y2); */
        
/*         sphere->x = x1; */
/*         sphere->y = y1; */
/*         if ( !test_sphere_safety(container, sphere) && state == 0  ){ */
/*             return NULL; */
/*         }else if (test_sphere_safety(container, sphere) && state == 1){ */
/*             state = 0; */
/*             printf("state = 0 left\n"); */
/*         } */
/*         sphere->z = sphere->z - delta_z;         */
/*     } */
    
    
/*     state = 0; */
/*     // initial position */
/*     sphere->z = z_max; */
/*     sphere->x = x_temp; */
/*     sphere->y = y_temp; */
/*     // test initial position */
/*     if (!test_sphere_safety(container, sphere)){ */
/*         printf("Initial position not safe!!!\n"); */
/*         state = 1; */
/*     } */
/*     sphere->z = sphere->z - delta_z; */
/*     // test right position     */
/*     while (sphere->z >= z_min){ */
/*         solve_x_y(sphere, first_sphere, second_sphere, sphere->z, &x1, &y1, &x2, &y2); */
        
/*         sphere->x = x2; */
/*         sphere->y = y2; */
/*         if ( !test_sphere_safety(container, sphere) && state == 0  ){ */
/*             return NULL; */
/*         }else if (test_sphere_safety(container, sphere) && state == 1){ */
/*             state = 0; */
/*             printf("state = 0 right\n"); */
/*         } */
/*         sphere->z = sphere->z - delta_z;         */
/*     } */

/*     printf("#########################################################################\n"); */
/*     printf("Something wrong here!!!\n"); */
/*     printf("#########################################################################\n"); */
/*     printf("z_max = %.3lf\n", z_max); */
/*     printf("z_min = %.3lf\n", z_min); */
/*     exit(1); */
/*     return NULL; */
/* } */



/* double compute_z_min(struct Sphere *s, struct Sphere *s1, struct Sphere *s2) */
/* { */
/*     double x1 = s1->x; */
/*     double y1 = s1->y; */
/*     double z1 = s1->z; */
/*     double r1 = s1->radius; */
    
/*     double x2 = s2->x; */
/*     double y2 = s2->y; */
/*     double z2 = s2->z; */
/*     double r2 = s2->radius;     */
    
/*     /\* double x = s->x; *\/ */
/*     /\* double y = s->y; *\/ */
/*     /\* double z = s->z; *\/ */
/*     double r = s->radius; */

/*     // First we transform everything into O1 frame */
/*     x2 = x2 - x1; */
/*     y2 = y2 - y1; */
/*     z2 = z2 - z1; */

/*     double dist = x2*x2+y2*y2+z2*z2; */
/*     dist = sqrt(dist); */
/*     if ( dist -(r1+r2+2.0*r) > 0.0001 ) */
/*         printf("Cannot do that!!!\n"); */

/*     double d = -z2/x2; */
/*     double e = (x2*x2 + y2*y2 + z2*z2 + (r+r1)*(r+r1) - (r+r2)*(r+r2))/(2.0*x2); */
/*     double p2 = y2/x2; */
/*     double a = p2*p2 + 1 + d*d; */
/*     double b = 2.0*d*e; */
/*     double c = e*e - (p2*p2 + 1)*(r+r1)*(r+r1); */
/*     // return the value in earth frame */
/*     if (b*b - 4.0*a*c < 0) */
/*         printf("NAN\n"); */

/*     return (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a) + z1;; */
/* } */

/* double compute_z_max(struct Sphere *s, struct Sphere *s1, struct Sphere *s2) */
/* { */
/*     double x1 = s1->x; */
/*     double y1 = s1->y; */
/*     double z1 = s1->z; */
/*     double r1 = s1->radius; */
    
/*     double x2 = s2->x; */
/*     double y2 = s2->y; */
/*     double z2 = s2->z; */
/*     double r2 = s2->radius;     */
    
/*     /\* double x = s->x; *\/ */
/*     /\* double y = s->y; *\/ */
/*     /\* double z = s->z; *\/ */
/*     double r = s->radius; */

/*     // First we transform everything into O1 frame */
/*     x2 = x2 - x1; */
/*     y2 = y2 - y1; */
/*     z2 = z2 - z1; */

/*     double dist = x2*x2+y2*y2+z2*z2; */
/*     dist = sqrt(dist); */
/*     if ( dist -(r1+r2+2.0*r) > 0.0001 ) */
/*         printf("Cannot do that!!!\n"); */

/*     double d = -z2/x2; */
/*     double e = (x2*x2 + y2*y2 + z2*z2 + (r+r1)*(r+r1) - (r+r2)*(r+r2))/(2.0*x2); */
/*     double p2 = y2/x2; */
/*     double a = p2*p2 + 1 + d*d; */
/*     double b = 2.0*d*e; */
/*     double c = e*e - (p2*p2 + 1)*(r+r1)*(r+r1); */
/*     // return the value in earth frame */
/*     if (b*b - 4.0*a*c < 0) */
/*         printf("NAN\n"); */

/*     return (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a) + z1;; */
/* } */


/* // */
/* //  Sphere s is contacting with both sphere s1 and s2. This function computes (x1, y1) and (x2, y2) */
/* //  given z. */
/* // */

/* void solve_x_y(struct Sphere* s, struct Sphere* s1, struct Sphere* s2, double z, double *x_1, double *y_1, double *x_2, double *y_2) */
/* { */
    
/*     double x1 = s1->x; */
/*     double y1 = s1->y; */
/*     double z1 = s1->z; */
/*     double r1 = s1->radius; */
    
/*     double x2 = s2->x; */
/*     double y2 = s2->y; */
/*     double z2 = s2->z; */
/*     double r2 = s2->radius;     */
    
/*     /\* double x = s->x; *\/ */
/*     /\* double y = s->y; *\/ */
/*     /\* double z = s->z; *\/ */
/*     double r = s->radius; */

/*     // First we transform everything into O1 frame */
/*     x2 = x2 - x1; */
/*     y2 = y2 - y1; */
/*     z2 = z2 - z1; */
/*     z = z - z1; */

/*     double p1 = (x2*x2 + y2*y2 + (r+r1)*(r+r1) - z*z - (r+r2)*(r+r2) + (z2-z)*(z2-z))/(2.0*x2); */
/*     double p2 = y2/x2; */
/*     double a = p2*p2+1; */
/*     double b = -2.0*p1*p2; */
/*     double c = p1*p1 - (r+r1)*(r+r1) + z*z; */

/*     *y_1 = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a); */
/*     *y_2 = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a); */
    
/*     *x_1 = p1 - p2*(*y_1); */
/*     *x_2 = p1 - p2*(*y_2); */

/*     // Then we transform back to the earth frame */
/*     *x_1 = *x_1 + x1; */
/*     *x_2 = *x_2 + x2; */

/*     *y_1 = *y_1 + y1; */
/*     *y_2 = *y_2 + y2; */
/* } */

