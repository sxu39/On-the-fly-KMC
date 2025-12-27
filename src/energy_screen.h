#ifndef ENERGY_SCREENING_H
#define ENERGY_SCREENING_H

#include <string>

using namespace std;

class EnergyScreen
{
    public:
        void energy_screen(string root_dir, string iniconfrlx_dir, string opt_folder_dir, string confcand_folder_dir, int conf_num, int max_iter, double upper_devi, double aveforce_limit);
        void extract_optenergy(string file_name, double& energy);
};
#endif
