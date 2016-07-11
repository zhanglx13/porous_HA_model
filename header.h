
#include <stdbool.h>  //bool

#define LENGTH   1500
#define HEIGHT   1500
#define WIDTH    1500

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
};

struct sphere_list
{
    struct Sphere *sphere;
    struct sphere_list *next;
};


double rand_ball(unsigned int, unsigned int, double, double);
void rand_init_position(double *, double *, struct Container, double);
void update_sphere(struct Sphere*, struct Container *);
void update_container(struct Sphere*, struct Container*);
bool collide_z(struct Sphere, struct Sphere);
void simulate(double, unsigned int, unsigned int, struct Container *);
struct Sphere * find_first_collision(struct Sphere *, struct Container *);
struct Sphere * find_second_collision(struct Shpere *);

/* helper functions */
bool test_overlap(struct Container);
bool overlap(struct Sphere, struct Sphere);
void logistic(struct Container);
void ball_info(struct Container, double, double, unsigned int*, unsigned int*);
double compute_z_collide(struct Sphere *, struct Sphere *);



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
 */
bool collide_z(struct Sphere s0, struct Sphere s1)
{
    double delta_x = s0.x - s1.x;
    double delta_y = s0.y - s1.y;
    return sqrt(delta_x * delta_x + delta_y * delta_y) < (s0.radius + s1.radius) ? true : false; 
}


bool overlap(struct Sphere s0, struct Sphere s1)
{
    double delta_x = s0.x - s1.x;
    double delta_y = s0.y - s1.y;
    double delta_z = s0.z - s1.z;
    return sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z) < (s0.radius + s1.radius) ? true : false;
}

bool test_overlap(struct Container container)
{
    struct sphere_list *first, *second;
    first = container.head;
    second = container.head;
    while(first != NULL){
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
    printf("Container's dim: (%.3lf, %.3lf, %.3lf)\n", container.length, container.width, container.height);
    printf("#Ball A: %u\n", numA);
    printf("#Ball B: %u\n", numB);
    printf("Ball positions:\n");
    struct sphere_list *ptr = container.head;
    if (ptr == NULL)
        printf("Empty container\n");
    while(ptr != NULL){
        printf("(%.3lf, %.3lf, %.3lf, %.3lf)\n", ptr->sphere->x, ptr->sphere->y, ptr->sphere->z, ptr->sphere->radius);
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
