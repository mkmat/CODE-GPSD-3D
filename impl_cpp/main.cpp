// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011
#include <iostream>
#include <vector>
#include "voro++/voro++.hh"
using namespace voro;
 
// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
 
// Set up the number of blocks that the container is divided into
const int n_x=1,n_y=1,n_z=1;
 
// Set the number of particles that are going to be randomly introduced
const int particles=20;
 
// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}
 
int main() {
 
         unsigned int i,j;
         int id,nx,ny,nz;
         double x,y,z;
         voronoicell c;
         std::vector<int> neigh,f_vert;
         std::vector<double> v;

         int num_faces;
         std::vector<double> f_areas;
         std::vector<int> f_orders;

         FILE *f = fopen("coords.csv", "w");
 
         // Create a container with the geometry given above, and make it
         // non-periodic in each of the three coordinates. Allocate space for
         // eight particles within each computational block
         container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                         false,false,false,particles);
 
         fprintf(f, "p_id,x,y,z\n");
         // Randomly add particles into the container
         for(i=0;i<particles;i++) {
                 x=x_min+rnd()*(x_max-x_min);
                 y=y_min+rnd()*(y_max-y_min);
                 z=z_min+rnd()*(z_max-z_min);
                 con.put(i,x,y,z);
                 //std::cout<<"\t x = "<<x<<"\t y = "<<"\t z = "<<z<<std::endl;
                 fprintf(f, "%d,%lf,%lf,%lf\n",i,x,y,z);

         }

         fclose(f);

         int idx;
         int f_sum = 0;

         FILE *verts = fopen("vertices.csv", "w");
         FILE *faces = fopen("faces.csv", "w");

         fprintf(verts, "p_id,v_id,x,y,z\n");
         fprintf(faces, "p_id,f_id,total,vertices\n");

         c_loop_all cl(con);
         int total_verts;

         if(cl.start()) do if(con.compute_cell(c,cl)) {
                 cl.pos(x,y,z);id=cl.pid();

                 //std::cout<<"id = "<<id<<"\t x = "<<x<<"\t y = "<<"\t z = "<<z<<std::endl;
 
                 // Gather information about the computed Voronoi cell
                 c.neighbors(neigh);
                 c.face_vertices(f_vert);
                 c.vertices(x,y,z,v);

                 c.face_areas(f_areas);
                 c.face_orders(f_orders);

                 //std::cout<<"vertices = "<<v.size()<<"\t actual = "<<v.size()/3<<std::endl;

                 f_sum = 0;
                 
                 for (int i = 0; i < f_orders.size(); i++){

                        total_verts = f_vert[f_sum];
                        fprintf(faces, "%d,%d,%d,[",id,i,total_verts);

                        for (int j = 0; j < total_verts; j++){

                                f_sum++;
                                
                                if (j != (f_orders[i] - 1))
                                        fprintf(faces, "%d.",   f_vert[f_sum]);
                                else
                                        fprintf(faces, "%d]\n", f_vert[f_sum]);
                        }

                        f_sum++;

                        //f_sum += f_orders[i];

                 }


                 idx = v.size()/3;

                 for (int i = 0; i < idx; i++){
                        fprintf(verts, "%d,%d,", id,i);

                        for (int j = 0; j < 3; j++){

                                if (j != 2)
                                        fprintf(verts, "%lf,",  v[3*i + j]);
                                else
                                        fprintf(verts, "%lf\n", v[3*i + j]);
                        }
                        
                 }




 

         } while (cl.inc());

        fclose(verts);
        fclose(faces);

 
}