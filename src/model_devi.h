#ifndef MODEL_DEVI_H
#define MODEL_DEVI_H

#include <vector>
#include <string>
#include "in_out.h"
#include "cell.h"
#include "deepmd/DeepPot.h"

using namespace std;

namespace KMC
{
    class ModelDevi
    {
	public:
            in_out io_obj;
            cell opt_cell;
        public:
	    ModelDevi(vector<string>& modelPaths);
            void cal_devi_m(string opt_folder_dir, int conf_num, double z_low, double z_high);
	    double cal_devi_s(string filename,double z_low, double z_high);
	private:
            vector<deepmd::DeepPot*> dp_models;
    };
}

#endif  // MODEL_DEVI_H
