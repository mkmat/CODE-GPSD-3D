#include <voronoi_faces.hh>

#ifndef VORONOI_PARTICLE_HH
#define VORONOI_PARTICLE_HH

namespace gpsd_3d
{

class voronoi_particle
{

public:

    void set_voronoi_face(std::vector<coords> vertices);
    void set_particle_coord(double x, double y, double z);
    void clear_faces();
    void print_particle();
    int num_faces();
    coords position;

private:
    std::vector<voronoi_faces> all_faces;

};

void voronoi_particle::clear_faces()
{
    all_faces.clear();
}

void voronoi_particle::set_voronoi_face(std::vector<coords> vertices)
{
    voronoi_faces temp;
    temp.set_voronoi_faces(vertices);
    all_faces.push_back(temp);
}

void voronoi_particle::set_particle_coord(double x, double y, double z)
{
    position.set_coords(x, y, z);
}

void voronoi_particle::print_particle()
{
    for (int i = 0; i < all_faces.size(); i++){
        std::cout<<"face "<<i<<std::endl;
        all_faces[i].print_face();
    }
}

int voronoi_particle::num_faces()
{
    return all_faces.size();
}


}


#endif