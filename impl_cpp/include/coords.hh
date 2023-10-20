#include <iostream>
#include <vector>

#ifndef COORDS_HH
#define COORDS_HH

namespace gpsd_3d{

class coords{

public:

    void set_coords(std::vector<double> v);
    void set_coords(double x_c, double y_c, double z_c);

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


}

#endif