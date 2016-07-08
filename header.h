
#include <stdbool.h>  //bool

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
void update_sphere(struct Sphere*, struct Container);
void update_container(struct Sphere, struct Container*);
bool collide_z(struct Sphere, struct Sphere);
