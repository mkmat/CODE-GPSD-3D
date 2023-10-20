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
    voronoi_faces(std::vector<std::vector<double>> vertices);
    ~voronoi_faces();
};

voronoi_faces::voronoi_faces(std::vector<std::vector<double>> vertices)
{

    num_vertices = vertices.size();

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

    

}

voronoi_faces::~voronoi_faces()
{
}



}

#endif