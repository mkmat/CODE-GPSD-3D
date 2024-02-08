// 10 nov 2023 mk@mat.ethz.ch

#include <iostream>
#include <vector>
#include <fstream>
#include "include-voro++.hh"
using namespace voro;
 
int main() {
 
         unsigned int i,j;
         int id,id_read;
         double x,y,z;
         voronoicell_neighbor c;
         std::vector<int> neigh,f_vert;
         std::vector<double> v;
         double x_min,x_max,y_min,y_max,z_min,z_max;
         double cvol;
         int voro_max_vertices=0;
         int voro_max_faces=0;
         int voro_max_vertices_for_face=0;
         int triangles=0;

         int num_faces;
         std::vector<double> f_areas;
         std::vector<int> f_orders;

         // read configuration
         id           = 0;
         std::ifstream myfile("config.txt",std::ios_base::in);
         if (myfile.is_open()) {
            while (myfile >> id_read) {
                myfile >> x;
                myfile >> y;
                myfile >> z;
                id += 1; 
            }
         }
         myfile.close();
         const int particles=id;

         // read box information
         std::ifstream boxfile("box.txt",std::ios_base::in);
         if (boxfile.is_open()) {
            boxfile >> x_min;
            boxfile >> x_max;
            boxfile >> y_min;
            boxfile >> y_max;
            boxfile >> z_min;
            boxfile >> z_max;
            cvol=(x_max-x_min)*(y_max-y_min)*(z_max-z_min);
         }
         boxfile.close();

         // number of blocks that the container is divided into
         double g = pow(cvol/particles,0.33333);
         int n_x = 1 + int((x_max-x_min)/g);
         int n_y = 1 + int((y_max-y_min)/g);
         int n_z = 1 + int((z_max-z_min)/g);

         // Create a container with the geometry given above, and make it
         // periodic in each of the three coordinates. Allocate space for
         // eight particles within each computational block
         container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                         true,true,true,particles);  

         // add particles to container
         std::ifstream copyfile("config.txt",std::ios_base::in);
         if (copyfile.is_open()) {
            while (copyfile >> id) {
                copyfile >> x;
                copyfile >> y;
                copyfile >> z;
                con.put(id,x,y,z);
            }
         }
         copyfile.close();

         int idx;

         c_loop_all cl(con);
         if(cl.start()) do if(con.compute_cell(c,cl)) {
                 cl.pos(x,y,z);id=cl.pid();

                 // Gather information about the computed Voronoi cell
                 c.neighbors(neigh);
                 c.face_vertices(f_vert);
                 c.vertices(x,y,z,v);
                 c.face_orders(f_orders);

                 idx = v.size()/3;
                 if (idx > voro_max_vertices) { voro_max_vertices = idx; }; 

                 int using_faces=0;
                 for (i=0, j=0; i < f_orders.size(); i++){
                    if (id>neigh[i]) {
                        using_faces += 1; 
                    }
                    j += f_vert[j]+1;
                 }
                 if (using_faces > voro_max_faces) { voro_max_faces = using_faces; }; 

                 for (i=0, j=0; i < f_orders.size(); i++){
                    if (id>neigh[i]) { 
                        triangles = triangles + f_vert[j];
                        if (f_vert[j] > voro_max_vertices_for_face) { voro_max_vertices_for_face = f_vert[j]; }; 
                    }
                    j += f_vert[j]+1; 
                 }

         } while (cl.inc());

         printf("%d %d %d %d %d\n",particles,voro_max_vertices,voro_max_faces,voro_max_vertices_for_face,triangles);

}
