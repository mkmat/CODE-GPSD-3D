#include <simulation_box.hh>
#include <chrono>
#include <ctime>

int main(int argc, char *argv[])
{
    std::chrono::time_point<std::chrono::system_clock> start, end;

    gpsd_3d::simulation_box test(argc, argv);
    test.calculate_gpsd();


    return 0;
}