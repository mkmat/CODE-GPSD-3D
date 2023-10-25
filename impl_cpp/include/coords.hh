#include <iostream>
#include <vector>
#include<cmath>

#ifndef COORDS_HH
#define COORDS_HH

namespace gpsd_3d{

class coords{

public:

    void set_coords(std::vector<double> v);
    void set_coords(double x_c, double y_c, double z_c);
    std::vector<double> convert_to_vector();
    void print_coords();
    double return_distance_sq(coords cx);
    double return_distance(coords cx);

    double x;
    double y;
    double z;

};

void coords::set_coords(std::vector<double> v)
{
    x = v[0];
    y = v[1];
    z = v[2];
}

void coords::set_coords(double x_c, double y_c, double z_c)
{
    x = x_c;
    y = y_c;
    z = z_c;
}

std::vector<double> coords::convert_to_vector()
{
    std::vector<double> temp;
    temp.resize(3);

    /*temp.push_back(x);
    temp.push_back(y);
    temp.push_back(z);*/

    temp[0] = x;
    temp[1] = y;
    temp[2] = z;

    return temp;
}

void coords::print_coords()
{
    std::cout<<x<<","<<y<<","<<z<<std::endl;
}

double coords::return_distance_sq(coords cx)
{
    return ((this->x-cx.x)*(this->x-cx.x) + (this->y-cx.y)*(this->y-cx.y) + (this->z-cx.z)*(this->z-cx.z));
}

double coords::return_distance(coords cx)
{
    return std::sqrt(return_distance_sq(cx));
}

}

#endif