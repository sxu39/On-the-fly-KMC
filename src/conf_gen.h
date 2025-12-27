#ifndef CONF_GEN_H
#define CONF_GEN_H

#include "cell.h"
#include "in_out.h"
#include "rand_seed.h"
#include <string> 

using namespace std;

namespace KMC
{
    class ConfGen
    {
	public:
	    in_out io_obj;
	    cell opt_cell;
	    RandSeed randseed;
        public:
	    void genconf(string root_dir, string gen_folder_dir,string input_file, double rcut, double rcut_fix, double radius, double z_fix, int n_iter, int num_yes, int conf_num);
    };
}

#endif  // CONF_GEN_H
