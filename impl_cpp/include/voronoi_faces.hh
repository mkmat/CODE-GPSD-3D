#include <triangle.hh>

#ifndef VORONOI_FACES_HH
#define VORONOI_FACES_HH

namespace gpsd_3d{

class voronoi_faces
{

public:
    void set_voronoi_faces(std::vector<coords> vertices);
    void print_face();

private:
    std::vector<triangle> face_triangles;
    coords geometrical_centre;
    int num_vertices;
};

void voronoi_faces::set_voronoi_faces(std::vector<coords> vertices)
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
            temp_triangle_vertices.push_back(vertices[i]);
            temp_triangle_vertices.push_back(vertices[i+1]);
            temp_triangle_vertices.push_back(geometrical_centre);
            temp_triangle.set_vertices(temp_triangle_vertices);
            face_triangles.push_back(temp_triangle);

        }

        temp_triangle_vertices.clear();
        temp_triangle_vertices.push_back(vertices[num_vertices-1]);
        temp_triangle_vertices.push_back(vertices[0]);
        temp_triangle_vertices.push_back(geometrical_centre);
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

}

void voronoi_faces::print_face()
{
    for (int i = 0; i < face_triangles.size(); i++){
        std::cout<<"triangle "<<i<<std::endl;
        face_triangles[i].print_triangle();
    }
}

}

#endif