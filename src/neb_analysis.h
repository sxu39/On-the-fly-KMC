#ifndef NEB_ANALYSIS_H
#define NEB_ANALYSIS_H

#include "rand_seed.h"
#include "copy.h"

using namespace std;
class NEBAna
{
    public:
        RandSeed randseed;
        Copy copy;
    public:
        void neb_analysis(string root_dir,string iniconfrlx_dir,string opt_folder_dir,string neb_folder_dir,int n_iter,int image_num,double temp,double lower_barrier,double upper_barrier);
};
#endif // NEB_PROFILE_H

