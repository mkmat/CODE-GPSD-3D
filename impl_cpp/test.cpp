#include <simulation_box.hh>

int main(int agrc, char *argv[])
{
    std::cout<<argv[0]<<"\t"<<argv[1]<<std::endl;
    gpsd_3d::simulation_box test(argv[1]);
    test.calculate_gpsd();
    //test.test_setup(1);

    /*gpsd_3d::coords a;
    gpsd_3d::coords b;

    a.set_coords(5.,3., 9.);
    b.set_coords(1.,1.,1.);
    b = a;


    b.print_coords();
    a.print_coords();*/

    return 0;


}