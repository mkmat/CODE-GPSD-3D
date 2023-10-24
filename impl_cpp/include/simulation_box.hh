#include<voronoi_particle.hh>
#include<voro++/voro++.hh>
#include<cstring>
#include<sstream>
#include<fstream>
#include<cmath>

#ifndef SIMULATION_BOX_HH
#define SIMULATION_BOX_HH

namespace gpsd_3d{

class simulation_box{

public:

    std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);    
    simulation_box(char *filename);
    void print_coords();
    void test_setup(int n);
    int  return_position_in_grid(int *pos); 

private:

    double ro = 1.;
    double rc = 0.;
    double rp = 0.;
    const int dim = 3;
    int num_particles;
    double r_max;
    double delta_x;
    int *nx;
    double xlo;
    double xhi;
    double ylo;
    double yhi;
    double zlo;
    double zhi;
    double *L;
    int    *nx;
    int    *L_eff;
    double *delta_x;
    double *inv_deltax;

    std::vector<voronoi_particle> all_particles;
    std::vector<coords> all_coords;
    std::vector<int> neighbour_ids;
    std::vector<std::vector<int>> grid;
    std::vector<int> neighbour_list;

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

    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,1,1,1,true,true,true,num_particles);

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

    r_max = 0.;

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

        //std::cout<<"faces = "<<c.number_of_faces()<<"\t"<<all_particles[id].num_faces()<<std::endl;
        if (c.max_radius_squared() > r_max)
            r_max = c.max_radius_squared();

    } while (cl.inc());

    r_max      = sqrt(r_max);

    L = (double*)malloc(sizeof(double) * dim);
    L[0] = xhi - xlo;
    L[1] = yhi - ylo;
    L[2] = zhi - zlo;

    delta_x    = (double*)malloc(sizeof(double) * dim);
    inv_deltax = (double*)malloc(sizeof(double) * dim);
    nx         = (int*)malloc(sizeof(int) * dim);

    for (int axis = 0; axis < dim; axis++){
        delta_x[axis]    = r_max;
        inv_deltax[axis] = 1./r_max;
        nx[axis]         = (int)(L[axis] * inv_deltax[axis]);
        delta_x[axis]    = L[axis]/(1. * nx[axis]);
        inv_deltax[axis] = 1. * delta_x[axis];
    }

    L_eff = (int*)malloc(sizeof(int)*dim);

    for (int axis = 0; axis < dim; axis++)
        L_eff[axis] = 1;

    for (int axis = 0; axis < dim; axis++)
    {
        for (int itr = axis+1; itr < dim; itr++)
            L_eff[axis] *= nx[itr];
    }

    int L_total = 1;

    for (int axis = 0; axis < dim; axis++)
        L_total *= nx[axis];

    grid.resize(L_total);

    int *position_in_grid;
    int *neigh_position_in_grid;

    position_in_grid       = (int*)malloc(sizeof(int) * dim);
    neigh_position_in_grid = (int*)malloc(sizeof(int) * dim);

    coords cx;
    int temp_index;

    for (int i = 0; i < num_particles; i++){

        cx = all_particles[i].position;

        position_in_grid[0] = (int)(cx.x * inv_deltax[0]);
        position_in_grid[1] = (int)(cx.y * inv_deltax[1]);
        position_in_grid[2] = (int)(cx.z * inv_deltax[2]);

        for (int ii = -1; ii <= 1; ii++){
            for (int jj = -1; jj <= 1; jj++){
                for (int kk = -1; kk <= 1; kk++){

                    neigh_position_in_grid[0] = position_in_grid[0] + ii;
                    neigh_position_in_grid[1] = position_in_grid[1] + jj;
                    neigh_position_in_grid[2] = position_in_grid[2] + kk;

                    for (int axis = 0; axis < dim; axis++){
                        temp_index = neigh_position_in_grid[axis];
                        neigh_position_in_grid[axis] += nx[axis] * ((temp_index >= nx[axis]) - (temp_index < 0));  
                    }

                    grid[return_position_in_grid(neigh_position_in_grid)].push_back(i);

                }
            }
        }


    }

    free(position_in_grid);
    free(neigh_position_in_grid);
    

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

int simulation_box::return_position_in_grid(int *pos)
{
    int counter = 0;

    for (int axis = 0; axis < dim; axis++)
        counter += pos[axis] * L_eff[axis];
}


}

#endif