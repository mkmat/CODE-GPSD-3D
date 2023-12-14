#include <iostream>
#include <vector>
#include<cmath>

#ifndef COORDS_HH
#define COORDS_HH

namespace gpsd_3d{

class coords{

public:

    coords(double cx, double cy, double cz);
    void set_coords(std::vector<double> v);
    void set_coords(double x_c, double y_c, double z_c);
    std::vector<double> convert_to_vector();
    void print_coords();
    double return_distance_sq(coords cx);
    double return_distance(coords cx);
    double return_norm();
    double return_norm_sq();
    void   normalize();
    void   calculate_norm();
    friend coords operator+(const coords &c1, const coords &c2);
    friend coords operator-(const coords &c1, const coords &c2);
    friend double operator*(const coords &c1, const coords &c2);
    friend coords operator*(const double &d1, const coords &c2);

    double x;
    double y;
    double z;
    double norm_sq;
    double norm;
    double norm_inv;

};

coords::coords(double cx=0., double cy=0., double cz=0.)
{
    x = cx; y = cy; z = cz;
}

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

double coords::return_norm_sq()
{
    norm = x*x + y*y + z*z;
    return norm;
}

double coords::return_norm()
{
    norm = std::sqrt(return_norm_sq());
    return norm;
}

void coords::calculate_norm()
{
    norm = std::sqrt(return_norm_sq());
}

void coords::normalize()
{
    norm_inv = 1./return_norm();

    x *= norm_inv;
    y *= norm_inv;
    z *= norm_inv;

}

coords operator+(const coords &c1, const coords &c2)
{
    return coords(c1.x+c2.x, c1.y+c2.y, c1.z+c2.z);
}

coords operator-(const coords &c1, const coords &c2)
{
    return coords(c1.x-c2.x, c1.y-c2.y, c1.z-c2.z);
}

double operator*(const coords &c1, const coords &c2)
{
    return (c1.x*c2.x + c1.y*c2.y + c1.z*c2.z);
}

coords operator*(const double &d1, const coords &c2)
{
    return coords(d1 * c2.x, d1 * c2.y, d1 * c2.z);
}


}

#endif