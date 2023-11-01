#include <triangle.hh>

#ifndef VORONOI_FACES_HH
#define VORONOI_FACES_HH

namespace gpsd_3d{

class voronoi_faces
{

public:
    void set_face_vertices(std::vector<coords> vertices);
    void set_normal(double nx, double ny, double nz);
    void print_face();
    double return_max_radius_for_face(coords p, coords vx, double rsc, coords &lpes_c);
    void set_new_coordinate_system(coords p);
    void print_n_coords();

private:
    std::vector<triangle> face_triangles;
    coords geometrical_centre;
    int num_vertices;
    std::vector<double> candidate_r_max;
    int num_triangles;
    double face_r_max;
    int    argmax;
    double triangle_r_max;
    coords material_projection;
    coords material_position;
    coords n_f;
    coords n_p;
    coords n_q;
    //coords test_n;
};

void voronoi_faces::set_face_vertices(std::vector<coords> vertices)
{

    num_vertices = vertices.size();
    face_triangles.clear();

    double gx = 0.;
    double gy = 0.;
    double gz = 0.;

    for (int i = 0; i < num_vertices; i++){

        gx += vertices[i].x;
        gy += vertices[i].y;
        gz += vertices[i].z;

    }

    gx /= (1. * num_vertices);
    gy /= (1. * num_vertices);
    gz /= (1. * num_vertices);

    geometrical_centre.set_coords(gx, gy, gz);

    triangle temp_triangle;
    std::vector<coords> temp_triangle_vertices;

    if (num_vertices > 3){

        for (int i = 0; i < (num_vertices-1); i++){

            temp_triangle_vertices.clear();
            temp_triangle_vertices.push_back(geometrical_centre);
            temp_triangle_vertices.push_back(vertices[i]);
            temp_triangle_vertices.push_back(vertices[i+1]);
            temp_triangle.set_vertices(temp_triangle_vertices);
            face_triangles.push_back(temp_triangle);

        }

        temp_triangle_vertices.clear();
        temp_triangle_vertices.push_back(geometrical_centre);
        temp_triangle_vertices.push_back(vertices[num_vertices-1]);
        temp_triangle_vertices.push_back(vertices[0]);
        temp_triangle.set_vertices(temp_triangle_vertices);
        face_triangles.push_back(temp_triangle);

    }

    else {
        temp_triangle_vertices.clear();            
        temp_triangle_vertices.push_back(vertices[0]);
        temp_triangle_vertices.push_back(vertices[1]);
        temp_triangle_vertices.push_back(vertices[2]);
        temp_triangle.set_vertices(temp_triangle_vertices);
        face_triangles.push_back(temp_triangle);
    }

    num_triangles = face_triangles.size();
}

void voronoi_faces::print_face()
{
    for (int i = 0; i < face_triangles.size(); i++){
        std::cout<<"triangle "<<i<<std::endl;
        face_triangles[i].print_triangle();
    }
}

void voronoi_faces::set_new_coordinate_system(coords p)
{
    n_p = p - (p*n_f) * n_f;
    n_p.normalize();

    n_q.set_coords((n_f.y*n_p.z)-(n_f.z*n_p.y), (n_f.z*n_p.x)-(n_f.x*n_p.z), (n_f.x*n_p.y)-(n_f.y*n_p.x));
    n_q.normalize();

    material_projection = (n_f * geometrical_centre) * n_f;

    //std::cout<<"dots = "<<face_triangles[0].vA * n_f<<"\t"<<face_triangles[0].vB * n_f<<"\t"<<face_triangles[0].vC * n_f<<std::endl;
    //std::cout<<"norms = "<<n_f*n_f<<"\t"<<n_p*n_p<<"\t"<<n_q*n_q<<std::endl; 
    //std::cout<<"dots  = "<<n_f*n_p<<"\t"<<n_p*n_q<<"\t"<<n_q*n_f<<std::endl; 
}

double voronoi_faces::return_max_radius_for_face(coords p, coords vx, double rs, coords &lpes_c)
{
    set_new_coordinate_system(p);
    face_r_max = 0.;

    for (int i = 0; i < num_triangles; i++){
        face_triangles[i].transform_all_coordinates(p, material_projection, n_f, n_p, n_q);
        triangle_r_max = face_triangles[i].return_max_distance_for_triangle(p, vx, rs, lpes_c);
        
        if (triangle_r_max > face_r_max)
            face_r_max = triangle_r_max;
    }

    return face_r_max;
}

void voronoi_faces::set_normal(double nx, double ny, double nz)
{
    n_f.set_coords(nx, ny, nz);
    n_f.normalize();

    //n_f.print_coords();
    //test_n.print_coords();
}

void voronoi_faces::print_n_coords()
{
    n_f.print_coords();
}

}

#endif