#include "conf_gen.h"
#include <iostream>
#include <string> 
#include <iomanip>
#include <algorithm>

void KMC::ConfGen::genconf(string root_dir, string gen_folder_dir,string input_file, double rcut, double rcut_fix,double radius, double z_fix, int n_iter, int num_yes, int conf_num)
{
    opt_cell = io_obj.read_stru(input_file);
    string head = io_obj.header;
    prec_t scale_factor = io_obj.scaling;
    vector<string> elem = io_obj.elements;
    int natoms = opt_cell.get_atom_num();
    const vector<prec_t> &mass = opt_cell.get_mass();
    vector<int> atom_sel(natoms);
    vector<index_t> atype = opt_cell.get_atype();
    array<double,9> box = opt_cell.get_box();
    vector<prec_t> ini_coords = opt_cell.positions.get_coords();
    vector<double> x_old,y_old,z_old,x_new,y_new,z_new,x_ini,y_ini,z_ini;
    x_ini.resize(natoms), y_ini.resize(natoms), z_ini.resize(natoms);
    x_old.resize(natoms), y_old.resize(natoms), z_old.resize(natoms);
    x_new.resize(natoms), y_new.resize(natoms), z_new.resize(natoms);
    for (int i = 0; i < natoms; i++)
    {
        x_old[i] = x_ini[i] = ini_coords[3*i];
        y_old[i] = y_ini[i] = ini_coords[3*i+1];
        z_old[i] = z_ini[i] = ini_coords[3*i+2];
    }
    double lx = box[0], ly = box[4], lz = box[8];
    auto maxElement = max_element(z_ini.begin(), z_ini.end());
    double max_z = *maxElement;
    double z_low,z_high;
    string upperpos_file = "upper_pos.log";
    string upperpos_path = root_dir + "/" + upperpos_file;
    ofstream ofs_pos(upperpos_path, ios::app);
    if (num_yes == 0)
    {
        z_low = z_fix;
        z_high = max_z - 3.0;     
	ofs_pos << n_iter - 1 << " " << z_low << " " << z_high << endl;
    }
    else if (num_yes > 0)
    {
	ifstream ifs_pos(upperpos_path, ios::in);
	int iter_id;
	double value1,value2;
	vector<double> higher_value;
	higher_value.clear();
	for (int i = 1; i <= n_iter - 1; i++)
	{
            ifs_pos >> iter_id >> value1 >> value2;
	    higher_value.push_back(value2);
	}
	string status_file = "status.log";
	string status_path = root_dir + "/" + status_file;
	ifstream ifs_status(status_path, ios::in);
	string line;
        string lastLine;
        while (getline(ifs_status, line))
        {
            lastLine = line;
        }
        istringstream iss(lastLine);
        int iternum;
	string str_yes;
	iss >> iternum >> value1;
        z_low = higher_value[iternum-1];
        z_high = max_z - 3.0;
	ofs_pos << n_iter << " " << z_low << " " << z_high << endl;
	ifs_pos.close();
	ifs_status.close();
    }
    string randlog_file = "random.log";
    string randlog_path = root_dir + "/" + randlog_file;
    ofstream ofs(randlog_path,ios::app);
    int n_count;
    int min = 0;
    int max = natoms - 1;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(min, max);
    int rand_num = dis(gen);
    while (z_ini[rand_num] < z_low || z_ini[rand_num]> z_high)
    {
        rand_num = dis(gen);
    }
    ofs << rand_num << endl;
    ofs.close();
    for (int i = 0; i < natoms; i++)
    {
        atom_sel[i] = 0;
        double dx = x_old[i] - x_old[rand_num];
        double dy = y_old[i] - y_old[rand_num];
        double dz = z_old[i] - z_old[rand_num];
        dx = dx - round(dx/lx)*lx;
        dy = dy - round(dy/ly)*ly;
        dz = dz - round(dz/lz)*lz;
        double dr = sqrt(dx*dx + dy*dy + dz*dz);
        if (dr < rcut && z_old[i] > z_fix)
        {
            atom_sel[i] = 1;
        }
    }
    for (int n_conf = 1; n_conf <= conf_num; n_conf++)
    {
        for (int i = 0; i < natoms; i++)
        {
            if (atom_sel[i] == 1)
            {
                bool flag = true;
                while(flag)
                {
                    double random_x = randseed.get_rand_seed();
                    double random_y = randseed.get_rand_seed();
                    double random_z = randseed.get_rand_seed();
                    x_new[i] = x_old[i] + (random_x - 0.5)*radius;
                    y_new[i] = y_old[i] + (random_y - 0.5)*radius;
                    z_new[i] = z_old[i] + (random_z - 0.5)*radius;
                    if (x_new[i] < 0.0) x_new[i] = x_new[i] + lx;
                    if (y_new[i] < 0.0) y_new[i] = y_new[i] + ly;
                    if (z_new[i] < 0.0) z_new[i] = z_new[i] + lz;
                    if (x_new[i] > lx) x_new[i] = x_new[i] - lx;
                    if (y_new[i] > ly) y_new[i] = y_new[i] - ly;
                    if (z_new[i] > lz) z_new[i] = z_new[i] - lz;
                    n_count = 0;
                    for (int j = 0; j < natoms; j++)
                    {
                        if (atom_sel[j] == 1)
                        {
                            x_ini[j] = x_new[j];
                            y_ini[j] = y_new[j];
                            z_ini[j] = z_new[j];
                        }
                        double dx = x_new[i] - x_ini[j];
                        double dy = y_new[i] - y_ini[j];
                        double dz = z_new[i] - z_ini[j];
                        dx = dx - round(dx/lx)*lx;
                        dy = dy - round(dy/ly)*ly;
                        dz = dz - round(dz/lz)*lz;
                        double dr = sqrt(dx*dx + dy*dy + dz*dz);
                        if (i!=j && dr >= 1.0)
                        {
                            n_count += 1;
                        }
                    }
		    if (n_count == (natoms - 1))
                    {
                        flag = false;
                    }
                }
            }
            else
            {
                x_new[i] = x_old[i];
                y_new[i] = y_old[i];
                z_new[i] = z_old[i];
            }
        }
	string output_file = gen_folder_dir + "/" + "trial_conf_" + to_string(n_conf);
        ofstream ofs(output_file);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs << head << "\n";
        ofs << setprecision(2) << scale_factor << "\n";
        ofs << setprecision(6);
        for (int i = 0 ; i < 3 ; i++)
            ofs << box[i*3] << " " << box[i*3+1] << " " << box[i*3+2] << "\n";
        for (auto ele : elem)
            ofs << ele << " ";
	ofs << "\n";
	unordered_map<index_t, size_t> ele_num;
	for (int i = 0 ; i < elem.size() ; i++)
            ele_num.insert(pair<index_t, size_t>(i, 0));
        for (auto type : atype)
            ele_num[type] += 1;
        for (int i = 0 ; i < elem.size() ; i++)
            ofs << ele_num[i] << " ";
        ofs << "\n";
	ofs << "Selective Dynamics\n";
	ofs << "Cartesian\n";
	for (int i = 0; i < natoms; i++)
        {
	    double dx = x_new[i] - ini_coords[3*rand_num];
            double dy = y_new[i] - ini_coords[3*rand_num+1];
            double dz = z_new[i] - ini_coords[3*rand_num+2];
            dx = dx - round(dx/lx)*lx;
            dy = dy - round(dy/ly)*ly;
            dz = dz - round(dz/lz)*lz;
            double dr = sqrt(dx*dx + dy*dy + dz*dz);
            if (dr <= rcut_fix && z_new[i] >= z_fix)
                ofs << setw(10)<< fixed << x_new[i] << setw(10)<< fixed << y_new[i] << setw(10)<< fixed << z_new[i] << " T" << " T" << " T" << endl;
	    else
 	        ofs << setw(10)<< fixed << x_new[i] << setw(10)<< fixed << y_new[i] << setw(10)<< fixed << z_new[i] << " F" << " F" << " F" << endl;
        }
	ofs.close();
    }
}
