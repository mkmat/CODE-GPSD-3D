#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    triangle(std::vector<std::vector<double>> vertices);
    triangle();
    bool check_viability();
    double calculate_max_distance();

private:

    coords vA;
    coords vB;
    coords vC;
    coords centroid;

};

triangle::triangle(std::vector<std::vector<double>> vertices){

    vA.set_coords(vertices[0]);
    vB.set_coords(vertices[1]);
    vC.set_coords(vertices[2]);

}

}

#endif