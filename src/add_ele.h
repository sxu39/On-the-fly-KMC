#ifndef ADD_ELE_H
#define ADD_ELE_H

#include <string>
#include "fileio.h"
#include "para.h"
#include "copy.h"
#include "cell.h"
#include "in_out.h"
#include "folder.h"
#include "verlet.h" 

using namespace std;

namespace KMC
{
    class AddEle
    {
        public:
            FileIO fileio;
            Para para;
            Copy copy;
            in_out io_obj;
            cell opt_cell_ini;
            cell opt_cell_ele;
	    Folder folder;
	    Verlet verlet;
        public:
            void addele(string root_dir,string iniconf_dir, string iniconfrlx_dir, int n_iter, int cont_num, int num_yes, int num_step,int dump_freq, int max_event, double temp); 
            void write_poscar(string iniconfrlx_path, string eleconf_path, string newconf_path, int n_iter);
    };
}
#endif
