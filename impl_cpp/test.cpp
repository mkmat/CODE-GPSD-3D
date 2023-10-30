#include <simulation_box.hh>

int main(int agrc, char *argv[])
{
    bool filename_b = false;
    bool outfile_b  = false;
    
    bool xlo_b = false;
    bool xhi_b = false;
    bool ylo_b = false;
    bool yhi_b = false;
    bool zlo_b = false;
    bool zhi_b = false;
    
    bool ro_b = false;
    bool rc_b = false;
    bool rp_b = false;
    bool q_b  = false;
    bool more_b = false;
    bool info_b = false;
    bool quiet_b = false;
    bool clean_b = false;

    char *in_filename;
    char *out_filename;
    
    double xlo;
    double xhi;
    double ylo;
    double yhi;
    double zlo;
    double zhi;
    
    double ro;
    double rp;
    double rc;

    

    


    std::vector<std::string> results;
    gpsd_3d::simulation_box test(argv[1]);


    return 0;
}