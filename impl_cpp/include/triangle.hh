#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    void set_vertices(std::vector<coords> vertices);
    void print_triangle();
    void get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p);
    double return_max_distance_for_triangle(coords p, coords vx, double rs);
    
    std::vector<double> distances;
    double r_max;
    double R1;
    double R2;
    double condition;

private:

    coords vA;
    coords vB;
    coords vC;
    coords centroid;

};

void triangle::set_vertices(std::vector<coords> vertices)
{
    vA.set_coords(vertices[0].x, vertices[0].y, vertices[0].z);
    vB.set_coords(vertices[1].x, vertices[1].y, vertices[1].z);
    vC.set_coords(vertices[2].x, vertices[2].y, vertices[2].z);

    double cx = 0.;
    double cy = 0.;
    double cz = 0.;

    cx = (vA.x + vB.x + vC.x)/3.;
    cy = (vA.y + vB.y + vC.y)/3.;
    cz = (vA.z + vB.z + vC.z)/3.;

    centroid.set_coords(cx, cy, cz);
}

void triangle::print_triangle()
{
    vA.print_coords();
    vB.print_coords();
    vC.print_coords();
}

void triangle::get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p)
{
    R1 = vx.return_distance(candidate_p) - rs;
    R2 = p.return_distance(candidate_p);
    std::cout<<"R1 = "<<R1<<"\t R2 = "<<R2<<std::endl;
    condition = 1.*(R1 > r_max)*(R1 >= R2);
    r_max = R1*condition + r_max*(1. - condition);
}

double triangle::return_max_distance_for_triangle(coords p, coords vx, double rs)
{
    r_max = 0.;

    get_candidate_rmax(p, vx, rs, vA);
    get_candidate_rmax(p, vx, rs, vB);
    get_candidate_rmax(p, vx, rs, vC);

    return r_max;
}

}

#endif