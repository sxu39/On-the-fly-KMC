#include "format.h"
#include <iostream>
#include <string> 
#include <iomanip>

void KMC::Format::pos2dump(string poscar_path, string dump_path)
{
    opt_cell = io_obj.read_stru(poscar_path);
    int natoms = opt_cell.get_atom_num();
    vector<index_t> atype = opt_cell.get_atype();
    array<double,9> box = opt_cell.get_box();
    vector<prec_t> coords = opt_cell.positions.get_coords();
    
    ofstream ofs(dump_path,ios::out);
    ofs << "ITEM: TIMESTEP\n";
    ofs << "0\n";
    ofs << "ITEM: NUMBER OF ATOMS\n";
    ofs << natoms << "\n";
    ofs << "ITEM: BOX BOUNDS pp pp pp\n";
    ofs << fixed << "0.000000" << " " << box[0] << endl;
    ofs << fixed << "0.000000" << " " << box[4] << endl;
    ofs << fixed << "0.000000" << " " << box[8] << endl;
    ofs << "ITEM: ATOMS id type x y z\n";
    for (int i = 0 ; i < natoms; i++)
    {
        ofs << fixed << i+1 << " " << atype[i]+1 << " " << coords[3*i]<< " " << coords[3*i+1] << " " << coords[3*i+2] << endl;
    }
}
