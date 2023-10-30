#include <coords.hh>

#ifndef TRIANGLE_HH
#define TRIANGLE_HH

namespace gpsd_3d{

class triangle{

public:

    void set_vertices(std::vector<coords> vertices);
    void print_triangle();
    void get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p, coords &lpes_c);
    void transform_all_coordinates(coords p, coords mp, coords n_f, coords n_p, coords n_q);
    coords transform_one_coordinate(coords cx, coords n_f, coords n_p, coords n_q);
    double return_max_distance_for_triangle(coords p, coords vx, double rs, coords &lpes_c);
    
    
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
    coords X_modified;
    coords P_modified;

    double a;
    double b;
    double c;
    double d;
    double e;

    double h;
    double f;
    double inv_C;

    double t_minus;
    double t_plus;
    double u_minus;
    double u_plus;
    double discriminant;
    double denominator;
    double numerator;

    coords analytical_r_max;
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

void triangle::get_candidate_rmax(coords p, coords vx, double rs, coords candidate_p, coords &lpes_c)
{
    R1 = vx.return_distance(candidate_p)-rs;
    R2 = p.return_distance(candidate_p);
    condition = 1.*(R1 > r_max)*(R1 >= R2);
    r_max  = R1*condition + r_max*(1. - condition);
    lpes_c = (condition * candidate_p) + ((1.-condition)*lpes_c);  
}

double triangle::return_max_distance_for_triangle(coords p, coords vx, double rs, coords &lpes_c)
{
    r_max = 0.;

    get_candidate_rmax(p, vx, rs, vA, lpes_c);
    get_candidate_rmax(p, vx, rs, vB, lpes_c);
    get_candidate_rmax(p, vx, rs, vC, lpes_c);

    inv_C = 1./vC_modified.z; 

    h = vA_modified.z*inv_C;
    f = vB_modified.z*inv_C;

    a = (2.*(vA_modified.y-(h*vC_modified.y))*P_modified.y) + (rs*rs) + (X_modified.x*X_modified.x) - (P_modified*P_modified);
    b = 2.*(vB_modified.y-(f*vC_modified.y))*P_modified.y;
    c = 4.*rs*rs*((X_modified.x*X_modified.x) + ((vA_modified.y-(h*vC_modified.y))*(vA_modified.y-(h*vC_modified.y))));
    d = 4.*rs*rs*((vA_modified.y-(h*vC_modified.y))*(vB_modified.y-(f*vC_modified.y)));
    e = 4.*rs*rs*((vB_modified - (f*vC_modified))*(vB_modified - (f*vC_modified)));

    discriminant = (a*a*e)-(2.*a*b*d)+(d*d)+(b*b*c)-(c*e);

    if (discriminant >= 0.)
    {
        numerator = (d - (a*b));
        discriminant = std::sqrt(discriminant);
        denominator = 1./((b*b) - e); 

        t_minus = (numerator - discriminant)*denominator;
        u_minus = -1.*(h+(f*t_minus));

        if ((t_minus >= 0.) && (u_minus > 0.) && ((t_minus+u_minus) <= 1.)){
            analytical_r_max = vA_modified + (t_minus * vB_modified) + (u_minus * vC_modified);
            R1 = vx.return_distance(analytical_r_max)-rs;
            r_max = R1*(R1 > r_max) + r_max*(r_max >= R1);
        }

        t_plus = (numerator + discriminant)*denominator;
        u_plus = -1.*(h+(f*t_plus));

        if ((t_plus >= 0.) && (u_plus > 0.) && ((t_plus+u_plus) <= 1.)){
            analytical_r_max = vA_modified + (t_plus * vB_modified) + (u_plus * vC_modified);
            R1 = vx.return_distance(analytical_r_max)-rs;
            r_max = R1*(R1 > r_max) + r_max*(r_max >= R1);
        }
    }

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
    P_modified     = transform_one_coordinate(p, n_f, n_p, n_q);
    origin         = transform_one_coordinate(mp, n_f, n_p, n_q);

    vA_modified       = vA_modified - origin;
    vB_modified       = vB_modified - origin;
    vC_modified       = vC_modified - origin;
    P_modified        = P_modified - origin;
    X_modified        = X_modified - origin;
    origin            = origin - origin; 

    vB_modified = vB_modified - vA_modified;
    vC_modified = vC_modified - vA_modified;   
}



}

#endif