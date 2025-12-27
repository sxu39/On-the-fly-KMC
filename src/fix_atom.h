#ifndef FIX_ATOM_H
#define FIX_ATOM_H

#include "cell.h"
#include "in_out.h"
#include <string>

using namespace std;

namespace KMC
{
    class FixAtom
    {
	public:
            in_out io_obj;
            cell opt_cell_ini;
	    cell opt_cell_fin;
        public:
	    void fix_atom(string initial_pos, string final_pos, string initial_pos_new, string final_pos_new, double maximum_dis);
	    void fix_image(string initial_pos, string final_pos, string initial_pos_new, string final_pos_new);
	    void fix_height(string initial_pos, string initial_pos_new, double z_fix);
    };
}

#endif  // FIX_ATOM_H
