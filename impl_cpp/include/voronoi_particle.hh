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
    coords position;
    double return_max_lpes_radius(coords p, double rs);

private:
    std::vector<voronoi_faces> all_faces;
    double particle_r_max;
    double face_r_max;
    int num_faces;

};

void voronoi_particle::clear_faces()
{
    all_faces.clear();
    num_faces = 0;
}

void voronoi_particle::set_voronoi_face(std::vector<coords> vertices)
{
    voronoi_faces temp;
    temp.set_voronoi_faces(vertices);
    all_faces.push_back(temp);
    num_faces += 1;
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


double voronoi_particle::return_max_lpes_radius(coords p, double rs)
{
    particle_r_max = 0.;

    for (int i = 0; i < num_faces; i++){
        face_r_max = all_faces[i].return_max_radius_for_face(p, position, rs);
        
        if (face_r_max > particle_r_max)
            particle_r_max = face_r_max;
    }

    return particle_r_max;

}


}


#endif