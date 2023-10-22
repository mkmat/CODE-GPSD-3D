#include <simulation_box.hh>

int main(int agrc, char *argv[])
{
    std::cout<<argv[0]<<"\t"<<argv[1]<<std::endl;
    gpsd_3d::simulation_box test(argv[1]);

}