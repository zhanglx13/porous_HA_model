
#include <stdbool.h>  //bool

#define LENGTH   1500
#define HEIGHT   1500
#define WIDTH    1500

#define NBALLS  300
#define RA       5
#define RB       50

#define NA      50
#define NB      100

#define RAND0(bound)        \
    (double)rand()*(double)(bound) / (double)RAND_MAX

struct Container
{
    double length;  // x coordinate
    double width;   // y coordinate
    double height;  // z coordinate
    struct sphere_list *head;   // First sphere in the container
};

struct Sphere
{
    double x;
    double y;
    double z;
    double radius;
    unsigned int id;
};

struct sphere_list
{
    struct Sphere *sphere;
    struct sphere_list *next;
};


#define HIT_BOTTOM     0
#define HIT_WALL_X     1
#define HIT_WALL_Y     2


double rand_ball(unsigned int, unsigned int, double, double);
void rand_init_position(double *, double *, struct Container, double);
void update_sphere(struct Sphere*, struct Container *);
void update_container(struct Sphere*, struct Container*);
bool collide_z(struct Sphere, struct Sphere);
void simulate(double, unsigned int, unsigned int, struct Container *);
struct Sphere * find_first_collision(struct Sphere *, struct Container *);
struct Sphere * find_second_collision(struct Sphere *, struct Sphere *, struct Container *, int *);
struct Sphere * find_third_collision();
double z_collide_two(struct Sphere *, struct Sphere *, double, double*, double*);

/* helper functions */
bool test_overlap(struct Container);
bool overlap(struct Sphere, struct Sphere);
void logistic(struct Container);
void ball_info(struct Container, double, double, unsigned int*, unsigned int*);
double compute_z_collide(struct Sphere *, struct Sphere *);
bool safe_container_boundary_sphere(struct Container *, struct Sphere *);
bool test_sphere_safety(struct Container *, struct Sphere *);



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

/*
 *   Test whether two spheres collide when projected in the x-y plane
 *   s1 should below s0
 */
bool collide_z(struct Sphere s0, struct Sphere s1)
{
    double delta_x = s0.x - s1.x;
    double delta_y = s0.y - s1.y;
    double dist = delta_x * delta_x + delta_y * delta_y - (s0.radius + s1.radius)*(s0.radius + s1.radius);
    return  ( dist < -0.00001 ) && (s1.z < s0.z)  ? true : false; 
}

/*
 *  Test if two balls collide
 */
bool overlap(struct Sphere s0, struct Sphere s1)
{
    double delta_x = s0.x - s1.x;
    double delta_y = s0.y - s1.y;
    double delta_z = s0.z - s1.z; 
    double dist = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
    double sum_radius = (s0.radius + s1.radius)*(s0.radius + s1.radius);
    double diff = dist - sum_radius;
    // Due to the existance of numerical error, the condition for overlap is released.
    if (diff < -0.0001){
        // printf("Ball %u (%.3lf, %.3lf, %.3lf) collides with ball %u (%.3lf, %.3lf, %.3lf)\n", s0.id, s0.x, s0.y, s0.z, s1.id, s1.x, s1.y, s1.z );
        // printf("dist = %.3lf\n", dist);
        // printf("sum_radius = %.3lf\n", sum_radius);
        // //exit(1);
        return true;
    }
    else{
        return false;
    }
    return  diff < - 0.0001 ? true : false;
}


/*
 *  Test if any two balls collide in the container
 */
bool test_overlap(struct Container container)
{
    struct sphere_list *first, *second;
    first = container.head;
    second = container.head;
    while(first != NULL){
        // first test container boudary safety
        if (!safe_container_boundary_sphere(&container, first->sphere))
            return true;
        // second test if two balls collide
        second = first->next;
        while(second != NULL){
            if (overlap(*(first->sphere), *(second->sphere)))
                return true;
            second = second->next;
        }
        first = first->next;
    }
    return false;
}


/*
 *  Test if the position of the sphere is safe in the container
 *  (The sphere is not included in the container yet)
 */
bool test_sphere_safety(struct Container *container, struct Sphere *sphere)
{
    if (!safe_container_boundary_sphere(container, sphere))
        return false;
    struct sphere_list *slist = container->head;
    while (slist != NULL){
        if (overlap(*sphere, *(slist->sphere) ))
            return false;
        slist = slist->next;
    }
    return true;
}


void update_container(struct Sphere *sphere, struct Container *container)
{
    struct sphere_list *new_sphere = (struct sphere_list*)malloc(sizeof(struct sphere_list));
    new_sphere->sphere = sphere;
    new_sphere->next = container->head;
    container->head = new_sphere;
    if (container->head == NULL)
        printf("Empty container\n");
    if (new_sphere == NULL)
        printf("Empty new sphere\n");
}



void logistic(struct Container container)
{
    unsigned int numA, numB;
    ball_info(container, RA, RB, &numA, &numB);
    printf("##########################################\n");
    printf("Container's dim: (%.3lf, %.3lf, %.3lf)\n", container.length, container.width, container.height);
    printf("#Ball A: %u\n", numA);
    printf("#Ball B: %u\n", numB);
    printf("Ball positions:\n");
    struct sphere_list *ptr = container.head;
    if (ptr == NULL)
        printf("Empty container\n");
    while(ptr != NULL){
        printf("ball %u: (%.3lf, %.3lf, %.3lf, %.3lf)\n", 
               ptr->sphere->id, ptr->sphere->x, ptr->sphere->y, ptr->sphere->z, ptr->sphere->radius);
        ptr = ptr->next;
    }
}

void ball_info(struct Container container, double Ra, double Rb, unsigned int *numA, unsigned int *numB)
{
    *numA = 0;
    *numB = 0;
    struct sphere_list *ptr = container.head;
    if(ptr == NULL)
        printf("Empty ptr\n");
    while(ptr != NULL){
        if (Ra == ptr->sphere->radius)
            *numA += 1;
        else if( Rb == ptr->sphere->radius)
            *numB += 1;
        else{ 
            printf("Error in ball dim: %.3lf\n",ptr->sphere->radius);         
            exit(0);
        }
        ptr = ptr->next;
    }
}

/* sphere s1 is moving */
double compute_z_collide(struct Sphere *s0, struct Sphere *s1)
{
    return  sqrt( (s0->radius + s1->radius)*(s0->radius + s1->radius) // (r+r`)^2
                - (s0->x - s1->x)*(s0->x - s1->x)                  //(x-x`)^2
                - (s0->y - s1->y)*(s0->y - s1->y) )                //(y-y`)^2
            + s0->z;
}

bool safe_container_boundary_sphere(struct Container *container, struct Sphere *sphere)
{
    return (   (sphere->x + sphere->radius <= container->length) && (sphere->x >= sphere->radius)
            && (sphere->y + sphere->radius <= container->width)  && (sphere->y >= sphere->radius)
            && (sphere->z + sphere->radius <= container->height) && (sphere->z >= sphere->radius)) 
            ? true : false;
}


/*   Given the position of two spheres s0 and s1, and the third sphere's radius
 *   this function computes the largest z coord when the third sphere collides
 *   with both s0 and s1.
 */
double z_collide_two(struct Sphere *s0, struct Sphere *s1, double r, double *x, double *y)
{
    double x0 = s0->x;
    double y0 = s0->y;
    double z0 = s0->z;
    double r0 = s0->radius;
    double x1 = s1->x;
    double y1 = s1->y;
    double z1 = s1->z;
    double r1 = s1->radius;

    // This problem is modeled as the following problem
    // Given two points of a triangle on the plane, also
    // given the three sides of the triangle, can we 
    // compute the coord of the third point?
    /*
                     C
                     /\
                    /  \
                   /    \
                  /      \
               A +--------+ B
        
        Given: A(0, 0), B(xb, yb), 
               |AB| = |(xb,yb)|
               |AC| = r0 + r
               |BC| = r1 + r
        Solve: C(xc, yc)

    */
    // Compute coord of B
    double xb, yb;
    xb = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
    yb = z1-z0;
    // Define the sides
    double AB = sqrt(xb*xb + yb*yb);
    double AC = r0+r;
    double BC = r1+r;
    
    /* printf("xb = %.3lf\n", xb); */
    /* printf("yb = %.3lf\n", yb); */
    
    /* printf("AB = %.3lf\n", AB); */
    /* printf("AC = %.3lf\n", AC); */
    /* printf("BC = %.3lf\n", BC); */
    
    // We first test if the triangle is a line    
    double dist = AC+BC-AB;
    if ( dist*dist < 0.0000001 )
        return (BC*z0 + AC*z1)/(AB);
    // Then we test if the given points can form a triangle
    if ( dist < -0.00001 )
        return -1;
    // Then we try to solve the triangle by solving the following equations:
    //  xc^2 + yc^2 = AC^2
    //  (xb-xc)^2 + (yb-yc)^2 = BC^2
    double p1 = (xb*xb+yb*yb+AC*AC-BC*BC)/(2.0*xb);
    double p2 = yb/xb;
    double a = p2*p2+1;
    double b = -2.0*p1*p2;
    double c = p1*p1-AC*AC;
    double yc = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
    double xc = p1-p2*yc;

    /* printf("xc = %.3lf\n", xc); */
    /* printf("yc = %.3lf\n", yc); */
    
    /* double err1 = AC*AC - xc*xc - yc*yc; */
    /* double err2 = BC*BC - (xb-xc)*(xb-xc) - (yb-yc)*(yb-yc); */
    /* printf("err1 = %.3lf\n", err1); */
    /* printf("err2 = %.3lf\n", err2); */

    //
    // Now compute the coordinates of the sphere in the absolute space
    //
    double p = (x1-x0)*xc/xb;
    double q = (y1-y0)*xc/xb;
    *x = x0 + p;
    *y = y0 + q;


    return yc + z0;
}
