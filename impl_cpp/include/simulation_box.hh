#include<voronoi_particle.hh>
#include<voro++/voro++.hh>
#include<cstring>
#include<sstream>
#include<fstream>

#ifndef SIMULATION_BOX_HH
#define SIMULATION_BOX_HH

namespace gpsd_3d{

class simulation_box{

public:

    std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);    
    simulation_box(char *filename);
    void print_coords();
    void test_setup(int n);

private:

    double L  = 10.;
    double ro = 1.;
    double rc = 0.;
    double rp = 0.;
    int num_particles;

    std::vector<voronoi_particle> all_particles;
    std::vector<coords> all_coords;

};

std::vector<std::string> simulation_box::split_string_by_delimiter(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

simulation_box::simulation_box(char *filename)
{

    std::ifstream parser(filename, std::ifstream::in);
    std::string str;
    std::vector<std::string> results;
    coords temp_c;

    num_particles = 0;

    while(getline(parser,str)){

        results = split_string_by_delimiter(str, ',');
        temp_c.set_coords(stod(results[0]), stod(results[1]), stod(results[2]));
        //con.put(num_particles++, temp_c.x, temp_c.y, temp_c.z);
        num_particles++;
        all_coords.push_back(temp_c);

    }

    parser.close();

    voro::container con(0.,L,0.,L,0.,L,1,1,1,true,true,true,num_particles);

    for (int i = 0; i < num_particles; i++){
        con.put(i, all_coords[i].x, all_coords[i].y, all_coords[i].z);
    }

    all_particles.resize(num_particles);
    voro::voronoicell c;
    std::vector<int> neigh,f_vert;
    std::vector<double> v;

    int num_faces;
    std::vector<double> f_areas;
    std::vector<int> f_orders;

    voro::c_loop_all cl(con);
    int total_verts;
    voronoi_particle temp_particle;
    double temp_x;
    double temp_y;
    double temp_z;
    int id;
    std::vector<coords> all_vertices;
    coords temp_vertex;
    int num_vertices;
    voronoi_faces temp_face;
    int f_sum;
    std::vector<coords> temp_face_vertex_coords;
    int v_idx;

    if(cl.start()) do if(con.compute_cell(c,cl)) {

        cl.pos(temp_x,temp_y,temp_z);
        id=cl.pid();

        temp_particle.clear_faces();
        temp_particle.set_particle_coord(temp_x, temp_y, temp_z);

        c.face_vertices(f_vert);
        c.vertices(temp_x, temp_y, temp_z, v);

        all_vertices.clear();
        num_vertices = v.size()/3;

        for (int i = 0; i < num_vertices; i++){
            temp_vertex.set_coords(v[3*i], v[3*i+1], v[3*i+2]);
            all_vertices.push_back(temp_vertex);
        }

        f_sum = 0;

        while (f_sum < f_vert.size()){

            total_verts = f_vert[f_sum];
            temp_face_vertex_coords.clear();

            //std::cout<<"hmmm"<<total_verts<<std::endl;

            for (int j = 0; j < total_verts; j++){
                f_sum += 1;
                v_idx  = f_vert[f_sum];
                //std::cout<<"v index "<<v_idx<<"\t"<<all_vertices.size()<<std::endl;
                temp_face_vertex_coords.push_back(all_vertices[v_idx]);
            }

            //std::cout<<f_sum<<"\t"<<f_vert.size()<<std::endl;

            temp_particle.set_voronoi_face(temp_face_vertex_coords);
            f_sum += 1;

        }

        all_particles[id] = temp_particle;

    } while (cl.inc());

}

void simulation_box::print_coords()
{
    for (int i = 0; i < all_coords.size(); i++){
        all_coords[i].print_coords();
    }
}

void simulation_box::test_setup(int n)
{
    all_particles[n].print_particle();
}


}

#endif