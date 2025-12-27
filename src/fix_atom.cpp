#include "fix_atom.h"
#include <iostream>
#include <string> 
#include <iomanip>

void KMC::FixAtom::fix_atom(string initial_pos, string final_pos, string initial_pos_new, string final_pos_new, double maximum_dis)
{
    opt_cell_ini = io_obj.read_stru(initial_pos);
    int natoms = opt_cell_ini.get_atom_num();
    vector<prec_t> coords_ini = opt_cell_ini.positions.get_coords();
    vector<double> x_ini,y_ini,z_ini;
    x_ini.resize(natoms), y_ini.resize(natoms), z_ini.resize(natoms);
    opt_cell_fin = io_obj.read_stru(final_pos);
    vector<prec_t> coords_fin = opt_cell_fin.positions.get_coords();
    array<double,9> box = opt_cell_fin.get_box();
    double lx = box[0], ly = box[4], lz = box[8];
    vector<double> x_fin,y_fin,z_fin;
    x_fin.resize(natoms), y_fin.resize(natoms), z_fin.resize(natoms);
    ifstream if_ini(initial_pos,ios::in);
    ifstream if_fin(final_pos,ios::in);
    ofstream of_ini(initial_pos_new,ios::out);
    ofstream of_fin(final_pos_new,ios::out);

    string line;
    while (getline(if_ini, line)) 
    {
	of_ini << line << endl;
        if (line.find("Direct") != string::npos || line.find("Cartesian") != string::npos) 
	{
            break;
        }
    }
    while (getline(if_fin, line))
    {
        of_fin << line << endl;
        if (line.find("Direct") != string::npos || line.find("Cartesian") != string::npos)
        {
            break;
        }
    }
    if_ini.close();
    if_fin.close();
    for (int i = 0; i < natoms; i++)
    {
        x_ini[i] = coords_ini[3*i];
        y_ini[i] = coords_ini[3*i+1];
        z_ini[i] = coords_ini[3*i+2];
        x_fin[i] = coords_fin[3*i];
        y_fin[i] = coords_fin[3*i+1];
        z_fin[i] = coords_fin[3*i+2];

        if (x_ini[i] < 0.0) x_ini[i] = x_ini[i] + lx;
        if (y_ini[i] < 0.0) y_ini[i] = y_ini[i] + ly;
        if (z_ini[i] < 0.0) z_ini[i] = z_ini[i] + lz;
        if (x_ini[i] > lx) x_ini[i] = x_ini[i] - lx;
        if (y_ini[i] > ly) y_ini[i] = y_ini[i] - ly;
        if (z_ini[i] > lz) z_ini[i] = z_ini[i] - lz;


        if (x_fin[i] < 0.0) x_fin[i] = x_fin[i] + lx;
        if (y_fin[i] < 0.0) y_fin[i] = y_fin[i] + ly;
        if (z_fin[i] < 0.0) z_fin[i] = z_fin[i] + lz;
        if (x_fin[i] > lx) x_fin[i] = x_fin[i] - lx;
        if (y_fin[i] > ly) y_fin[i] = y_fin[i] - ly;
        if (z_fin[i] > lz) z_fin[i] = z_fin[i] - lz;

        double dx = x_fin[i] - x_ini[i];
        double dy = y_fin[i] - y_ini[i];
        double dz = z_fin[i] - z_ini[i];
        double dr = sqrt(dx*dx + dy*dy + dz*dz);
	if (dr > maximum_dis)
	{
	    of_ini << setw(15)<< fixed << x_ini[i] << setw(15)<< fixed << y_ini[i] << setw(15)<< fixed << z_ini[i] << " T" << " T" << " T" << endl;
	    of_fin << setw(15)<< fixed << x_fin[i] << setw(15)<< fixed << y_fin[i] << setw(15)<< fixed << z_fin[i] << " T" << " T" << " T" << endl;
	}          
	else       
	{          
	    of_ini << setw(15)<< fixed << x_ini[i] << setw(15)<< fixed << y_ini[i] << setw(15)<< fixed << z_ini[i] << " F" << " F" << " F" << endl;
	    of_fin << setw(15)<< fixed << x_fin[i] << setw(15)<< fixed << y_fin[i] << setw(15)<< fixed << z_fin[i] << " F" << " F" << " F" << endl;
	}              
    }
    of_ini.close();
    of_fin.close();
}

void KMC::FixAtom::fix_image(string initial_pos, string final_pos, string initial_pos_new, string final_pos_new)
{
    opt_cell_ini = io_obj.read_stru(initial_pos);
    opt_cell_fin = io_obj.read_stru(final_pos);
    int natoms = opt_cell_ini.get_atom_num();
    vector<prec_t> coords_ini = opt_cell_ini.positions.get_coords();
    vector<prec_t> coords_fin = opt_cell_fin.positions.get_coords();
    ifstream if_ini(initial_pos,ios::in);
    ifstream if_fin(final_pos,ios::in);
    ofstream of_ini(initial_pos_new,ios::out);
    ofstream of_fin(final_pos_new,ios::out);

    string line;
    while (getline(if_ini, line))
    {
        of_ini << line << endl;
        if (line.find("Direct") != string::npos || line.find("Cartesian") != string::npos)
        {
            break;
        }
    }
    while (getline(if_fin, line))
    {
        of_fin << line << endl;
        if (line.find("Direct") != string::npos || line.find("Cartesian") != string::npos)
        {
            break;
	}
    }
    if_ini.close();
    if_fin.close();
    for (int i = 0; i < natoms; i++)
    {
        of_ini << setw(15)<< fixed << coords_ini[3*i] << setw(15)<< fixed << coords_ini[3*i+1] << setw(15)<< fixed << coords_ini[3*i+2] << " F" << " F" << " F" << endl;
	of_fin << setw(15)<< fixed << coords_fin[3*i] << setw(15)<< fixed << coords_fin[3*i+1] << setw(15)<< fixed << coords_fin[3*i+2] << " F" << " F" << " F" << endl;
    }
    of_ini.close();
    of_fin.close();
}

void KMC::FixAtom::fix_height(string initial_pos,string initial_pos_new,double z_fix)
{
    opt_cell_ini = io_obj.read_stru(initial_pos);
    int natoms = opt_cell_ini.get_atom_num();
    vector<prec_t> coords_ini = opt_cell_ini.positions.get_coords();
    ifstream if_ini(initial_pos,ios::in);
    ofstream of_ini(initial_pos_new,ios::out);

    string line;
    while (getline(if_ini, line))
    {
        of_ini << line << endl;
        if (line.find("Direct") != string::npos || line.find("Cartesian") != string::npos)
        {
            break;
        }
    }
    if_ini.close();
    for (int i = 0; i < natoms; i++)
    {
	double z_value = coords_ini[3*i+2];
	if (z_value < z_fix)
            of_ini << setw(15)<< fixed << coords_ini[3*i] << setw(15)<< fixed << coords_ini[3*i+1] << setw(15)<< fixed << coords_ini[3*i+2] << " F" << " F" << " F" << endl;
        else
	    of_ini << setw(15)<< fixed << coords_ini[3*i] << setw(15)<< fixed << coords_ini[3*i+1] << setw(15)<< fixed << coords_ini[3*i+2] << " T" << " T" << " T" << endl;
    }
    of_ini.close();
}
