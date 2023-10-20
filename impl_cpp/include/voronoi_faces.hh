#include <triangle.hh>

#ifndef VORONOI_FACES_HH
#define VORONOI_FACES_HH

namespace gpsd_3d{

class voronoi_faces
{
private:
    std::vector<triangle> face_triangles;
    coords geometrical_centre;

public:
    voronoi_faces(std::vector<std::vector<double>> vertices);
    ~voronoi_faces();
};

voronoi_faces::voronoi_faces(std::vector<std::vector<double>> vertices)
{

    double gx;
    double gy;
    double gz;

    

}

voronoi_faces::~voronoi_faces()
{
}



}

#endif