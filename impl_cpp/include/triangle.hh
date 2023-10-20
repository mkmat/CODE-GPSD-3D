#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    //triangle(std::vector<std::vector<double>> vertices);
    //triangle();
    void set_vertices(std::vector<std::vector<double>> vertices);
    bool check_viability();
    double calculate_max_distance();

private:

    coords vA;
    coords vB;
    coords vC;
    coords centroid;

};

void triangle::set_vertices(std::vector<std::vector<double>> vertices){

    vA.set_coords(vertices[0]);
    vB.set_coords(vertices[1]);
    vC.set_coords(vertices[2]);

    double cx = 0.;
    double cy = 0.;
    double cz = 0.;

    for (int i = 0; i < 3; i++){
        cx 
    }

}

}

#endif