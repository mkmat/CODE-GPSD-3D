#include <iostream>
#include <vector>

#ifndef COORDS_HH
#define COORDS_HH

namespace gpsd_3d{

class coords{

public:

    void set_coords(std::vector<double> v);

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


}

#endif