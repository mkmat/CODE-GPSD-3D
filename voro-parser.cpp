// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

// modified : 13 oct 2023 mk@mat.ethz.ch 
//          : 20 oct 2023 double-checked against command-line utility

// compile via: g++ -I /usr/local/include/voro++/ voro-parser.cpp /usr/local/lib/libvoro++.a -o voro-parser.ex

#include <iostream>
#include <vector>
#include <fstream>
#include "include-voro++.hh"
using namespace voro;
 
// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}
 
int main() {
 
         unsigned int i,j;
         int id,id_read;
         double x,y,z;
         voronoicell_neighbor c;
         std::vector<int> neigh,f_vert;
         std::vector<double> v;
         double x_min,x_max,y_min,y_max,z_min,z_max;
         double cvol;
         // int max_vertices,max_faces,max_vertices_for_face;
         // long int triangles

         int num_faces;
         std::vector<double> f_areas;
         std::vector<int> f_orders;
         // std::vector<double> f_normals;     // MK 

         // max_vertices = 0;
         // max_faces    = 0;
         // max_vertices_for_face = 0;
         // triangles    = 0;

         // determine number of particles
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
         printf("%d\n",particles);

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
                         true,true,true,particles);  // MK: set up container    MK false -> true

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

         // con.print_custom("%i %x %y %z %s %w %P %t","packing.custom");  // REFERENCE STATEMENT

         int idx;

         c_loop_all cl(con);
         if(cl.start()) do if(con.compute_cell(c,cl)) {
                 cl.pos(x,y,z);id=cl.pid();

                 // Gather information about the computed Voronoi cell
                 c.neighbors(neigh);
                 c.face_vertices(f_vert);
                 c.vertices(x,y,z,v);

                 // taken out: https://math.lbl.gov/voro++/examples/polygons/
                 //c.face_areas(f_areas);
                 c.face_orders(f_orders);
                 //c.normals(f_normals);

                 idx = v.size()/3;
                 printf("%lf %lf %lf ",x,y,z);   // (x,y,z) 
                 printf("%d ",idx);              // #vertices for id

                 // obtain MK_used number of faces
                 int MK_using_faces=0;
                 for (i=0, j=0; i < f_orders.size(); i++){
                    if (id>neigh[i]) {
                        MK_using_faces += 1; 
                        // if (f_vert[j]>max_vertices_for_face) { max_vertices_for_face=f_vert[j]; };
                    }
                    j += f_vert[j]+1;
                 }

                 // printf("%d ",f_orders.size());  // #faces for id
                 printf("%d ",MK_using_faces);  // #faces of id (the used ones only)

                 //if (idx>max_vertices) { max_vertices=idx; }; 
                 //if (f_orders.size()>max_faces) { max_faces = f_orders.size(); }; 
                 //if (MK_using_faces>max_faces) { max_faces = MK_using_faces; };

                 // ADD: list of the number of vertices for all #faces faces
                 // ADD: but take out those that are otherwise occurring twice
                 // the used faces are marked 
                 for (i=0, j=0; i < f_orders.size(); i++){
                    // ith face of particle id 
                    // f_vert[..] (0) number of faces f (1) face 1 (2) ... (f+1) number of faces (.) ...
                    // id corresponding corresponding to ith face is neigh[i]
                    if (id>neigh[i]) { 
                        printf("%d ",f_vert[j]); // #ids for faces
                    }
                    j += f_vert[j]+1; 
                 }
                 printf("\n");

                 // ADD: all #vertices vertex coordinates for vertex i
                 for (i=0; i < idx; i++){
                        printf("%lf %lf %lf ",v[3*i],v[3*i+1],v[3*i+2]); // vertex coordinates for vertex i
                 }

                 // ADD: all vertex ids for all faces, one face after the other
                 for (i=0, j=0; i < f_orders.size(); i++){
                        if (id>neigh[i]) {
                            // triangles += f_vert[j];
                            for (int k = 0; k < f_vert[j]; k++){
                                printf("%d ",1+f_vert[j+k+1]);  // vertex id for face 
                            }
                        }
                        j += f_vert[j]+1; 
                 }
                 printf("\n");    


         } while (cl.inc());

         // FILE *allocate = fopen("allocate.txt", "w");         // MK ADDED
         // fprintf(allocate, "%d %d %d %d",max_vertices,max_faces,max_vertices_for_face,triangles);
         // fclose(allocate);


}
