#include <simulation_box.hh>

int main(int agrc, char *argv[])
{
    std::cout<<argv[0]<<"\t"<<argv[1]<<std::endl;
    gpsd_3d::simulation_box test(argv[1]);
    test.calculate_gpsd();

    /*gpsd_3d::coords a(1,1,1);
    gpsd_3d::coords b(2,-1,-2);
    gpsd_3d::coords c;

    c = a-b;
    c.print_coords();

    double dot = 3;
    c = dot*b;
    //std::cout<<dot<<std::endl;
    c.print_coords();*/

    return 0;


}