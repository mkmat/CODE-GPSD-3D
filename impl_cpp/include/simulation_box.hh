#include<voronoi_particle.hh>
#include<voro++/voro++.hh>
#include<cstring>
#include<sstream>
#include<fstream>
#include<random>

#ifndef SIMULATION_BOX_HH
#define SIMULATION_BOX_HH

namespace gpsd_3d{

class simulation_box{

public:

    friend std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);    
    simulation_box(char *filename);
    void print_coords();
    void test_setup(int n);
    int  return_position_in_grid(int *pos); 
    void get_position_in_grid(coords cx);
    void calculate_gpsd();
    void get_LPES();
    bool check_probe_centre_viability();
    coords get_probe_centre_image(coords a, coords b);

private:

    double ro = 1.;
    double rc = 0.;
    double rp = 0.;
    double rs;
    double r_coated;
    const int dim = 3;
    int    num_shots = 1000;
    int    num_particles;
    double r_max;
    double r_max_squared;
    int    *nx;
    double xlo;
    double xhi;
    double ylo;
    double yhi;
    double zlo;
    double zhi;
    double *L;
    double *inv_L;
    int    *L_eff;
    double *delta_x;
    double *inv_deltax;
    int    *position_in_grid;    
    double periodic_distance_sq;
    double x_periodic;
    double y_periodic;
    double z_periodic;
    double diff;
    coords lpes_c;

    std::vector<voronoi_particle> all_particles;
    std::vector<coords> all_coords;
    std::vector<int> neighbour_ids;
    std::vector<std::vector<int>> grid;
    std::vector<int> neighbour_list;
    std::mt19937 generator;
    std::uniform_real_distribution<double> dis;
    std::vector<int> temp_neighbour_list;  

    coords lpes_centre;
    coords probe_centre;
    coords probe_centre_image;
    coords temp_periodic_coords;

};

std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter)
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

    xlo = 0.;
    xhi = 10.;
    ylo = 0.;
    yhi = 10.;
    zlo = 0.;
    zhi = 10.;

    rs = ro + rc + rp;
    r_coated = ro + rc;

    std::ifstream parser(filename, std::ifstream::in);
    std::string str;
    std::vector<std::string> results;
    coords temp_c;

    num_particles = 0;
    generator.seed(1729);

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
    int f_count;
    std::vector<coords> temp_face_vertex_coords;
    int v_idx;
    std::vector<double> face_normals;
    coords temp_r_max;

    r_max = 0.;

    if(cl.start()) do if(con.compute_cell(c,cl)) {

        cl.pos(temp_x,temp_y,temp_z);
        id=cl.pid();

        temp_particle.clear_faces();
        temp_particle.set_particle_coord(temp_x, temp_y, temp_z);

        c.face_vertices(f_vert);
        //c.vertices(temp_x, temp_y, temp_z, v);
        c.vertices(v);
        c.normals(face_normals);

        all_vertices.clear();
        num_vertices = v.size()/3;

        for (int i = 0; i < num_vertices; i++){
            temp_vertex.set_coords(v[3*i], v[3*i+1], v[3*i+2]);
            all_vertices.push_back(temp_vertex);
        }

        f_sum = 0;
        f_count = 0;

        while (f_sum < f_vert.size()){

            total_verts = f_vert[f_sum];
            temp_face_vertex_coords.clear();

            for (int j = 0; j < total_verts; j++){
                f_sum += 1;
                v_idx  = f_vert[f_sum];
                temp_face_vertex_coords.push_back(all_vertices[v_idx]);
            }

            temp_particle.set_voronoi_faces(temp_face_vertex_coords, face_normals[3*f_count], face_normals[3*f_count+1], face_normals[3*f_count+2]);
            f_sum   += 1;
            f_count += 1;

            /*temp_particle.set_face_normal(face_normals[3*f_count], face_normals[3*f_count+1], face_normals[3*f_count+2]);
            f_count += 1;*/

        }

        all_particles[id] = temp_particle;

        //std::cout<<"faces = "<<c.number_of_faces()<<"\t"<<all_particles[id].num_faces<<std::endl;
        for (int i = 0; i < num_vertices; i++){
            temp_r_max.set_coords(v[3*i], v[3*i+1], v[3*i+2]);

            if (temp_r_max.return_norm_sq() > r_max)
                r_max = temp_r_max.return_norm_sq();
        }


    } while (cl.inc());

    /*for (int i = 0; i < num_particles; i++)
    {
        std::cout<<"--------------------------\n";
        all_particles[i].print_normal();
        std::cout<<"--------------------------\n";
        exit(EXIT_FAILURE);
    }*/

    r_max      = std::sqrt(r_max);
    //std::cout<<"r_max = "<<r_max<<std::endl;
    //exit(EXIT_FAILURE);

    L = (double*)malloc(sizeof(double) * dim);
    L[0] = xhi - xlo;
    L[1] = yhi - ylo;
    L[2] = zhi - zlo;

    inv_L = (double*)malloc(sizeof(double) * dim);

    for (int axis = 0; axis < dim; axis++)
        inv_L[axis] = 1./L[axis];

    delta_x    = (double*)malloc(sizeof(double) * dim);
    inv_deltax = (double*)malloc(sizeof(double) * dim);
    nx         = (int*)malloc(sizeof(int) * dim);

    for (int axis = 0; axis < dim; axis++){
        delta_x[axis]    = r_max;
        inv_deltax[axis] = 1./r_max;
        nx[axis]         = (int)(L[axis] * inv_deltax[axis]);
        delta_x[axis]    = L[axis]/(1. * nx[axis]);
        inv_deltax[axis] = 1./delta_x[axis];
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

    //for (int axis = 0; axis < dim; axis++)
        //std::cout<<L_eff[axis]<<std::endl;

    //std::cout<<"L total = "<<L_total<<"grid size = "<<grid.size()<<std::endl;


    int *neigh_position_in_grid;

    position_in_grid       = (int*)malloc(sizeof(int) * dim);
    neigh_position_in_grid = (int*)malloc(sizeof(int) * dim);

    coords cx;
    int temp_index;

    for (int i = 0; i < num_particles; i++){

        cx = all_particles[i].position;
        get_position_in_grid(cx);

        //cx.print_coords();

        //for (int axis = 0; axis < dim; axis++)
            //std::cout<<"first = "<<inv_deltax[axis]<<std::endl;

        for (int ii = -1; ii <= 1; ii++){
            for (int jj = -1; jj <= 1; jj++){
                for (int kk = -1; kk <= 1; kk++){

                    neigh_position_in_grid[0] = position_in_grid[0] + ii;
                    neigh_position_in_grid[1] = position_in_grid[1] + jj;
                    neigh_position_in_grid[2] = position_in_grid[2] + kk;

                    for (int axis = 0; axis < dim; axis++){
                        temp_index = neigh_position_in_grid[axis];
                        neigh_position_in_grid[axis] += nx[axis] * ((temp_index < 0) - (temp_index >= nx[axis]));
                        //std::cout<<"temp_index = "<<temp_index<<" axis = "<<axis<<" index = "<<neigh_position_in_grid[axis]<<std::endl;
                    }

                    //std::cout<<"position in grid = "<<return_position_in_grid(neigh_position_in_grid)<<std::endl;
                    grid[return_position_in_grid(neigh_position_in_grid)].push_back(i);

                }
            }
        }


    }

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

    return counter;
}

void simulation_box::get_position_in_grid(coords cx)
{
    position_in_grid[0] = (int)(cx.x * inv_deltax[0]);
    position_in_grid[1] = (int)(cx.y * inv_deltax[1]);
    position_in_grid[2] = (int)(cx.z * inv_deltax[2]);
}

void simulation_box::calculate_gpsd()
{
    bool condition;
    int  num_count = 0;

    double particle_max;
    double lpes_max;
    int i_index=1729;
    FILE *f;
    f = fopen("r_max_details.csv", "w");
    fprintf(f, "id,px,py,pz,cx,cy,cz,lpes,type\n");
    std::string sol_type;


    while (num_count < 10000){

        probe_centre.set_coords(xlo+L[0]*dis(generator), ylo+L[1]*dis(generator), zlo+L[2]*dis(generator));
        //probe_centre.set_coords(7.25927,5.26002,3.81528);
        condition = check_probe_centre_viability();

        //probe_centre = all_particles[0].position;
        //condition = check_probe_centre_viability();

        if(condition){

            //probe_centre.print_coords();
            lpes_max = 0.;

            //std::cout<<"------------------------\n";

            for (int i = 0; i < temp_neighbour_list.size(); i++){

                probe_centre_image = get_probe_centre_image(probe_centre, all_particles[temp_neighbour_list[i]].position);
                //std::cout<<probe_centre_image.return_norm()<<"\t"<<probe_centre.return_distance(all_particles[temp_neighbour_list[i]].position)<<std::endl;
                //std::cout<<probe_centre_image.return_norm()<<"\t"<<r_max<<std::endl;

                if (probe_centre_image.return_norm() <= r_max){

                    //std::cout<<"heeereee"<<std::endl;
                    /*std::cout<<"-------------------------\n";
                    std::cout<<"probe centre       = "; probe_centre.print_coords();
                    std::cout<<"probe centre image = "; probe_centre.print_coords();
                    std::cout<<"material position  = "; all_particles[temp_neighbour_list[i]].position.print_coords();
                    std::cout<<"-------------------------\n";*/


                    particle_max = all_particles[temp_neighbour_list[i]].return_max_lpes_radius(probe_centre_image, rs, lpes_c, sol_type, lpes_max);
                    //particle_max = all_particles[0].return_max_lpes_radius(probe_centre, rs);
                    //std::cout<<"particle max = "<<particle_max<<std::endl;
                    if (particle_max > lpes_max){
                        lpes_max = particle_max;
                        i_index = i;
                    }

                }
            }

            //std::cout<<"r_max = "<<lpes_max<<" "<<i_index<<" "; probe_centre.print_coords();
            //std::cout<<lpes_max<<"\n";

            //probe_centre.print_coords();

            fprintf(f, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%s\n", num_count, probe_centre.x, probe_centre.y, probe_centre.z, lpes_c.x, lpes_c.y, lpes_c.z, lpes_max, sol_type.c_str());
            num_count++;

            //probe_centre.print_coords();

        }

    }

    fclose(f);
}

bool simulation_box::check_probe_centre_viability()
{
    temp_neighbour_list.clear();
    get_position_in_grid(probe_centre);
    temp_neighbour_list = grid[return_position_in_grid(position_in_grid)];

    for (int i = 0; i < temp_neighbour_list.size(); i++){
        probe_centre_image = get_probe_centre_image(probe_centre, all_particles[temp_neighbour_list[i]].position);
        if (probe_centre_image.return_norm() <= r_coated)
            return false;
    }

    return true;
}

coords simulation_box::get_probe_centre_image(coords a, coords b)
{
    diff = (a.x - b.x);
    x_periodic = diff - (L[0]*round(diff*inv_L[0]));

    diff = (a.y - b.y);
    y_periodic = diff - (L[1]*round(diff*inv_L[1]));

    diff = (a.z - b.z);
    z_periodic = diff - (L[2]*round(diff*inv_L[2]));

    temp_periodic_coords.set_coords(x_periodic, y_periodic, z_periodic);

    return temp_periodic_coords;
}


}

#endif