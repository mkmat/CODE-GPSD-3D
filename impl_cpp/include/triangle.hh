#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    //triangle(std::vector<std::vector<double>> vertices);
    //triangle();
    void set_vertices(std::vector<coords> vertices);
    void print_triangle();
    //bool check_viability();
    //double calculate_max_distance();

private:

    coords vA;
    coords vB;
    coords vC;
    coords centroid;

};

void triangle::set_vertices(std::vector<coords> vertices){

    /*vA.set_coords(vertices[0]);
    vB.set_coords(vertices[1]);
    vC.set_coords(vertices[2]);*/

    vA.set_coords(vertices[0].x, vertices[0].y, vertices[0].z);
    vB.set_coords(vertices[1].x, vertices[1].y, vertices[1].z);
    vC.set_coords(vertices[2].x, vertices[2].y, vertices[2].z);

    double cx = 0.;
    double cy = 0.;
    double cz = 0.;

    //for (int i = 0; i < 3; i++){
    cx = (vA.x + vB.x + vC.x)/3.;
    cy = (vA.y + vB.y + vC.y)/3.;
    cz = (vA.z + vB.z + vC.z)/3.;
    //}

    centroid.set_coords(cx, cy, cz);

}

void triangle::print_triangle()
{
    vA.print_coords();
    vB.print_coords();
    vC.print_coords();
}

}

#endif