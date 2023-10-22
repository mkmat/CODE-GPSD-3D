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

    while(getline(parser,str)){

        results = split_string_by_delimiter(str, ',');
        temp_c.set_coords(stod(results[0]), stod(results[1]), stod(results[2]));
        all_coords.push_back(temp_c);

    }

    parser.close();

    num_particles = all_coords.size();
    voro::voronoicell c;
    std::vector<int> neigh,f_vert;
    std::vector<double> v;

    int num_faces;
    std::vector<double> f_areas;
    std::vector<int> f_orders;

    //FILE *f = fopen("coords.csv", "w");

    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block
    voro::container con(0.,L,0.,L,0.,L,1,1,1,
                    true,true,true,num_particles);


    for (int i = 0; i < num_particles; i++){

    }

}

void simulation_box::print_coords()
{
    for (int i = 0; i < all_coords.size(); i++){
        all_coords[i].print_coords();
    }
}


}

#endif