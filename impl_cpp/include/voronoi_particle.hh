#include <voronoi_faces.hh>

#ifndef VORONOI_PARTICLE_HH
#define VORONOI_PARTICLE_HH

namespace gpsd_3d
{

class voronoi_particle
{

public:

    void set_voronoi_face(std::vector<std::vector<double>> vertices);
    void clear_faces();

private:
    coords position;
    std::vector<voronoi_faces> all_faces;

};

void voronoi_particle::clear_faces()
{
    all_faces.clear();
}

void voronoi_particle::set_voronoi_face(std::vector<std::vector<double>> vertices)
{
    voronoi_faces temp;
    temp.set_voronoi_faces(vertices);
    all_faces.push_back(temp);
}



}


#endif