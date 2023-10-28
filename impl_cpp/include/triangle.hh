#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    void set_vertices(std::vector<coords> vertices);
    void print_triangle();
    void get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p);
    void transform_all_coordinates(coords p, coords mp, coords n_f, coords n_p, coords n_q);
    coords transform_one_coordinate(coords cx, coords n_f, coords n_p, coords n_q);
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

    coords vA_modified;
    coords vB_modified;
    coords vC_modified;
    coords origin;
    coords material_position;
    coords probe_position;
};

void triangle::set_vertices(std::vector<coords> vertices)
{
    vA.set_coords(vertices[0].x, vertices[0].y, vertices[0].z);
    vB.set_coords(vertices[1].x, vertices[1].y, vertices[1].z);
    vC.set_coords(vertices[2].x, vertices[2].y, vertices[2].z);

}

void triangle::print_triangle()
{
    vA.print_coords();
    vB.print_coords();
    vC.print_coords();
}

void triangle::get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p)
{
    R1 = vx.return_distance(candidate_p);
    R2 = p.return_distance(candidate_p);
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

coords triangle::transform_one_coordinate(coords cx, coords n_f, coords n_p, coords n_q)
{
    return coords(cx*n_f, cx*n_p, cx*n_q);
}

void triangle::transform_all_coordinates(coords p, coords mp, coords n_f, coords n_p, coords n_q)
{
    vA_modified    = transform_one_coordinate(vA, n_f, n_p, n_q);
    vB_modified    = transform_one_coordinate(vB, n_f, n_p, n_q);
    vC_modified    = transform_one_coordinate(vC, n_f, n_p, n_q);
    probe_position = transform_one_coordinate(p, n_f, n_p, n_q);
    origin         = transform_one_coordinate(mp, n_f, n_p, n_q);

    vA_modified       = vA_modified - origin;
    vB_modified       = vB_modified - origin;
    vC_modified       = vC_modified - origin;
    probe_position    = probe_position - origin;
    material_position = material_position - origin;
    origin            = origin - origin;
    
}



}

#endif