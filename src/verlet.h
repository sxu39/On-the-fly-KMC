#ifndef VERLET_H
#define VERLET_H

#include <string>
#include "cell.h"
#include "in_out.h"
#include "common.h"

using namespace std;

namespace KMC
{
    class Verlet
    {
        public:
            in_out io_obj;
            cell opt_cell;
        public:
            void conf_md(string root_dir,string subiniconf_dir, int num_step, int dump_freq, double temp);
	    vector<double> gen_velocity(double temp, double mass);
	    vector<double> lj_force(vector<double> z_coord, double z_wall, int natoms);
    };
}
#endif
