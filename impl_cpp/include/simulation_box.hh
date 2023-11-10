#include<voronoi_particle.hh>
#include<voro++/voro++.hh>
#include<cstring>
#include<sstream>
#include<fstream>
#include<random>
#include <chrono>
#include <ctime>


#ifndef SIMULATION_BOX_HH
#define SIMULATION_BOX_HH

namespace gpsd_3d{

class simulation_box{

public:

    friend std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);    
    simulation_box(int argc, char *argv[]);
    void print_coords();
    void test_setup(int n);
    int  return_position_in_grid(int *pos); 
    void get_position_in_grid(coords cx);
    void calculate_gpsd();
    void print_info_file();
    void print_info_numbers_file();     
    bool check_probe_centre_viability();
    coords get_probe_centre_image(coords a, coords b);

private:

    double ro;
    double rc;
    double rp;
    double rs;
    double r_coated;
    const int dim = 3;
    int    num_shots;
    int    num_count = 0;  
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
    bool   more_b;
    bool   info;

    int    total_shots;
    double abs_lpes_min;
    double abs_lpes_max;
    int    total_num_triangles = -1;     // MK question: where is this defined or calculated?
    double pore_mean;
    double pore_std;

    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;

    double voro_time;
    double triangle_time;
    double mc_time;

    std::vector<voronoi_particle> all_particles;
    std::vector<coords> all_coords;
    std::vector<int> neighbour_ids;
    std::vector<std::vector<int>> grid;
    std::vector<int> neighbour_list;
    std::mt19937 generator;
    std::uniform_real_distribution<double> dis;
    std::vector<int> temp_neighbour_list; 


    std::string out_filename;
    std::string str_delimiter;   
    char mydelimiter = ',';     

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

simulation_box::simulation_box(int argc, char *argv[])
{
    bool xlo_b = false;
    bool xhi_b = false;
    bool ylo_b = false;
    bool yhi_b = false;
    bool zlo_b = false;
    bool zhi_b = false;
    bool filename_b = false;
    double q   = 10.0;         
    bool outfilename_b = false;
    bool boxname_b = false;    

    more_b = false;
    info   = false;
    
    std::string filename;
    std::string boxname;     

    ro = 1.;
    rc = 0.;
    rp = 0.;

    std::vector<std::string> results;

    for (int i = 1; i < argc; i++)
    {
        results.clear();
        results = gpsd_3d::split_string_by_delimiter(argv[i],'=');

        if (results[0] == "-in"){
            filename   = results[1];
            filename_b = true;
        }

        if (results[0] == "-box"){          
            boxname    = results[1];       
            boxname_b  = true;             
        }                                  

        if (results[0] == "-xlo"){
            xlo   = std::stod(results[1]);
            xlo_b = true;
        }

        if (results[0] == "-xhi"){
            xhi   = std::stod(results[1]);
            xhi_b = true;
        }

        if (results[0] == "-ylo"){
            ylo   = std::stod(results[1]);
            ylo_b = true;
        }

        if (results[0] == "-yhi"){
            yhi   = std::stod(results[1]);
            yhi_b = true;
        }

        if (results[0] == "-zlo"){
            zlo   = std::stod(results[1]);
            zlo_b = true;
        }

        if (results[0] == "-zhi"){
            zhi   = std::stod(results[1]);
            zhi_b = true;
        }

        if (results[0] == "-ro"){
            ro = std::stod(results[1]);
        }

        if (results[0] == "-rp"){
            rp = std::stod(results[1]);
        }

        if (results[0] == "-rc"){
            rc = std::stod(results[1]);
        }

        if (results[0] == "-q"){
            q  = std::stod(results[1]);  
        }

        if (results[0] == "-o"){
            out_filename = results[1];
            outfilename_b = true;
        }

        if (results[0] == "-more")
            more_b = true;

        if (results[0] == "-info")
            info = true;

        if (results[0] == "-d")                     
            str_delimiter = results[1];             
            mydelimiter = str_delimiter[1];         
        
    }

    std::cout << "delimiter [" << mydelimiter << "]\n";

    if (!filename_b)   
        std::cout<<"configuration file name -in=<filename> is missing"<<std::endl;

    if (boxname_b == true) {  
        std::ifstream myfile (boxname);
        if (myfile.is_open()) {
            std::string str;
            std::getline (myfile,str);
            results = split_string_by_delimiter(str, mydelimiter); 
            xlo = std::stod(results[0]);
            xhi = std::stod(results[1]);
            ylo = std::stod(results[2]);
            yhi = std::stod(results[3]);
            zlo = std::stod(results[4]);
            zhi = std::stod(results[5]);
            myfile.close();
        } else { 
            std::cout<<"unable to open boxfile " << boxname <<std::endl;
            exit(EXIT_FAILURE);
        };
    } else { 
        if (!xlo_b || !xhi_b || !ylo_b || !yhi_b || !zlo_b || !zhi_b) { 
            std::cout<<"one of the necessary arguments is missing"<<std::endl;
            exit(EXIT_FAILURE);
        }; 
    };

    rs = ro + rc + rp;
    r_coated = ro + rc;

    if (!outfilename_b)
        out_filename = "results.gpsd";
   
    std::ifstream parser(filename, std::ifstream::in);
    std::string str;
    results.clear();
    coords temp_c;

    num_particles = 0;
    generator.seed(1729);

    while(getline(parser,str)){
        results = split_string_by_delimiter(str, mydelimiter); 
        temp_c.set_coords(stod(results[1]), stod(results[2]), stod(results[3]));
        num_particles++;
        all_coords.push_back(temp_c);
    }
    parser.close();

    num_shots = q * num_particles; 


    start = std::chrono::high_resolution_clock::now();

    voro::container con(xlo,xhi,ylo,yhi,zlo,zhi,1,1,1,true,true,true,num_particles);

    for (int i = 0; i < num_particles; i++){
        con.put(i, all_coords[i].x, all_coords[i].y, all_coords[i].z);
    }

    end = std::chrono::high_resolution_clock::now();

    voro_time  = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    voro_time *= 1e-9; 

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

    start = std::chrono::high_resolution_clock::now();

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
        }

        all_particles[id] = temp_particle;

        for (int i = 0; i < num_vertices; i++){
            temp_r_max.set_coords(v[3*i], v[3*i+1], v[3*i+2]);

            if (temp_r_max.return_norm_sq() > r_max)
                r_max = temp_r_max.return_norm_sq();
        }


    } while (cl.inc());

    end = std::chrono::high_resolution_clock::now();

    triangle_time  = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    triangle_time *= 1e-9;     


    r_max      = 2. * std::sqrt(r_max);

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

    int *neigh_position_in_grid;

    position_in_grid       = (int*)malloc(sizeof(int) * dim);
    neigh_position_in_grid = (int*)malloc(sizeof(int) * dim);

    coords cx;
    int temp_index;

    for (int i = 0; i < num_particles; i++){

        cx = all_particles[i].position;
        get_position_in_grid(cx);

        for (int ii = -1; ii <= 1; ii++){
            for (int jj = -1; jj <= 1; jj++){
                for (int kk = -1; kk <= 1; kk++){

                    neigh_position_in_grid[0] = position_in_grid[0] + ii;
                    neigh_position_in_grid[1] = position_in_grid[1] + jj;
                    neigh_position_in_grid[2] = position_in_grid[2] + kk;

                    for (int axis = 0; axis < dim; axis++){
                        temp_index = neigh_position_in_grid[axis];
                        neigh_position_in_grid[axis] += nx[axis] * ((temp_index < 0) - (temp_index >= nx[axis]));
                    }

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
    double particle_max;
    double lpes_max;
    FILE *f;
    f = fopen(out_filename.c_str(), "w");

    if (more_b)
        fprintf(f, "id,px,py,pz,cx,cy,cz,r\n");

    std::string sol_type;

    total_shots  = 0;
    pore_mean    = 0.;
    pore_std     = 0.;
    abs_lpes_max = 0.;
    abs_lpes_min = L[0]*L[0]+L[1]*L[1]+L[2]*L[2];

    start = std::chrono::high_resolution_clock::now();


    while (total_shots < num_shots){ 

        probe_centre.set_coords(xlo+L[0]*dis(generator), ylo+L[1]*dis(generator), zlo+L[2]*dis(generator));
        condition = check_probe_centre_viability();

        if(condition){

            lpes_max = 0.;

            for (int i = 0; i < num_particles; i++){

                probe_centre_image = get_probe_centre_image(probe_centre, all_particles[i].position);
                particle_max = all_particles[i].return_max_lpes_radius(probe_centre_image, rs, lpes_c, sol_type, lpes_max);

                if (particle_max > lpes_max){
                    lpes_max = particle_max;
                }

            }

            if (more_b)
                fprintf(f, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", num_count, probe_centre.x, probe_centre.y, probe_centre.z, lpes_c.x, lpes_c.y, lpes_c.z, lpes_max);
            else
                fprintf(f, "%lf\n", lpes_max); 

            num_count++;
            pore_mean += lpes_max;
            pore_std  += (lpes_max * lpes_max);

            abs_lpes_max = (abs_lpes_max * (abs_lpes_max > lpes_max)) + (lpes_max * (lpes_max > abs_lpes_max));
            abs_lpes_min = (abs_lpes_min * (abs_lpes_min < lpes_max)) + (lpes_max * (lpes_max < abs_lpes_min));

        }

        total_shots += 1;

    }

    fclose(f);

    end      = std::chrono::high_resolution_clock::now();
    mc_time  = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    mc_time *= 1e-9;

    pore_mean = pore_mean/(1.*num_count); 
    pore_std  = pore_std/(1.*num_count);  
    pore_std  = pore_std - (pore_mean*pore_mean);

    if (info) {
        print_info_file();
        print_info_numbers_file(); 
    };

}

void simulation_box::print_info_file()
{
    std::string info_filename = "";
    std::vector<std::string> results;

    results = split_string_by_delimiter(out_filename, '.');

    for (int i = 0; i < (results.size()-1); i++)
        info_filename += results[i];

    info_filename += "-INFO"; 
    info_filename += "."+results[(results.size()-1)];

    FILE *f;
    f = fopen(info_filename.c_str(), "w");

    fprintf(f, "N=%d\n",num_particles);
    fprintf(f, "ro=%lf\n",ro);
    fprintf(f, "rp=%lf\n",rp);
    fprintf(f, "rc=%lf\n",rc);
    fprintf(f, "V=%lf\n", L[0] * L[1] * L[2]);                 
    fprintf(f, "triangles = %d\n", total_num_triangles);
    fprintf(f, "shots=%d\n", total_shots);
    fprintf(f, "radii=%d\n", num_count); 
    fprintf(f, "min_pore_radius=%lf\n", abs_lpes_min);
    fprintf(f, "max_pore_radius=%lf\n", abs_lpes_max);
    fprintf(f, "mean_pore_radius=%lf\n", pore_mean);    
    fprintf(f, "stderr_pore_radius=%lf\n", pore_std);
    fprintf(f, "phi_reff=%lf\n",1.0-float(num_count)/float(total_shots)); 
    fprintf(f, "cells=1\n");
    fprintf(f, "threads=1\n"); 
    fprintf(f, "voro_cpu_time=%lf\n", voro_time);    
    fprintf(f, "voro_real_time=%lf\n", voro_time);   
    fprintf(f, "triangle_setup_cpu_time=%lf\n", triangle_time); 
    fprintf(f, "triangle_setup_real_time=%lf\n", triangle_time);
    fprintf(f, "MonteCarlo_cpu_time=%lf\n", mc_time); 
    fprintf(f, "MonteCarlo_real_time=%lf\n", mc_time);

    fclose(f);

}

void simulation_box::print_info_numbers_file() 
{
    std::string info_filename = "";
    std::vector<std::string> results;

    results = split_string_by_delimiter(out_filename, '.');

    for (int i = 0; i < (results.size()-1); i++)
        info_filename += results[i];

    info_filename += "-INFO-numbers";   
    info_filename += "."+results[(results.size()-1)];

    FILE *f;
    f = fopen(info_filename.c_str(), "w");

    fprintf(f, "%d\n",num_particles);
    fprintf(f, "%lf\n",ro);
    fprintf(f, "%lf\n",rp);
    fprintf(f, "%lf\n",rc);
    fprintf(f, "%lf\n", L[0] * L[1] * L[2]);
    fprintf(f, "%d\n", total_num_triangles);
    fprintf(f, "%d\n", total_shots);
    fprintf(f, "%d\n", num_count);               
    fprintf(f, "%lf\n", abs_lpes_min);
    fprintf(f, "%lf\n", abs_lpes_max);
    fprintf(f, "%lf\n", pore_mean);
    fprintf(f, "%lf\n", pore_std);
    fprintf(f, "%lf\n", 1.0-float(num_count)/float(total_shots)); 
    fprintf(f, "1\n");  // number of neighbor cells
    fprintf(f, "1\n");  // number of threads used in parallel
    fprintf(f, "%lf\n", voro_time); // voro cpu time 
    fprintf(f, "%lf\n", voro_time); // voro real time
    fprintf(f, "%lf\n", triangle_time); // triangle cpu time
    fprintf(f, "%lf\n", triangle_time); // triangle real time
    fprintf(f, "%lf\n", mc_time);   // MC cpu time
    fprintf(f, "%lf\n", mc_time);   // MC real time
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
