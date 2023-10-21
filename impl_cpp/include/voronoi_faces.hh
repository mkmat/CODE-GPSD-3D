#include <triangle.hh>

#ifndef VORONOI_FACES_HH
#define VORONOI_FACES_HH

namespace gpsd_3d{

class voronoi_faces
{
private:
    std::vector<triangle> face_triangles;
    coords geometrical_centre;
    int num_vertices;

public:
    void set_voronoi_faces(std::vector<std::vector<double>> vertices);
};

void voronoi_faces::set_voronoi_faces(std::vector<std::vector<double>> vertices)
{

    num_vertices = vertices.size();
    face_triangles.clear();

    double gx = 0.;
    double gy = 0.;
    double gz = 0.;

    for (int i = 0; i < num_vertices; i++){

        gx += vertices[i][0];
        gy += vertices[i][1];
        gz += vertices[i][2];

    }

    gx /= (1. * num_vertices);

    geometrical_centre.set_coords(gx, gy, gz);
    std::vector<double> gc_vec = geometrical_centre.convert_to_vector();

    std::vector<std::vector<double>> temp_triangle_vertices;
    triangle temp_triangle;

    for (int i = 0; i < (num_vertices-1); i++){

        temp_triangle_vertices.clear();
        
        temp_triangle_vertices.push_back(vertices[i]);
        temp_triangle_vertices.push_back(vertices[i+1]);
        temp_triangle_vertices.push_back(gc_vec);

        temp_triangle.set_vertices(temp_triangle_vertices);
        face_triangles.push_back(temp_triangle);

    }

    temp_triangle_vertices.clear();
    temp_triangle_vertices.push_back(vertices[num_vertices-1]);
    temp_triangle_vertices.push_back(vertices[0]);
    temp_triangle_vertices.push_back(gc_vec);

    temp_triangle.set_vertices(temp_triangle_vertices);
    face_triangles.push_back(temp_triangle);

}

}

#endif