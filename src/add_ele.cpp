#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include "add_ele.h"

using namespace std;

void KMC::AddEle::addele(string root_dir, string iniconf_dir, string iniconfrlx_dir, int n_iter, int cont_num, int num_yes,int num_step,int dump_freq, int max_event, double temp)
{   
    int eventnum,total_eventnum,sum_eventnum=0;
    vector <int> event_num;
    string iter_file = "iter.log";
    string iter_path = root_dir + "/" + iter_file;
    ifstream ifs(iter_path,ios::in);
    string line,value;
    event_num.reserve(n_iter);
    for (int iter = 1; iter <= n_iter; iter++)
    {   
        for (int j = 0; j < 8; j++)
        {
            getline(ifs,line);
            if (line.find("Event_num") != string::npos)
            {
                stringstream ss(line);
                ss >> value >> eventnum;
                event_num.push_back(eventnum);
        	if (eventnum > 0)
		{	
		    if (num_yes == 0)
                    {
        	        sum_eventnum += 1;
		    }
                    else if (num_yes > 0)
                    {
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
                        iss >> iternum >> str_yes; 
			if (iter > iternum)
                        {
			    sum_eventnum += 1;
			} 
		    }
	        }
            }             
        }                 
    }                     
    cout << "sum event is " << sum_eventnum << endl; 
    ifs.close();          
    total_eventnum = 0;   
    for (int i = n_iter -  cont_num; i < n_iter; i++)
    {                     
        total_eventnum +=  event_num[i];
    }                     
    if (total_eventnum ==  0 || sum_eventnum >= max_event)
    {
        string iniconfrlx_file = "ini_conf_rlx";
        string iniconfrlx_path = iniconfrlx_dir + "/" + iniconfrlx_file;
        string subiniconf_folder = to_string(num_yes+2) + "Layer";
        string subiniconf_dir = iniconf_dir + "/" + subiniconf_folder;
        const char* subiniconf_folderpath = subiniconf_dir.c_str();
        if (access(subiniconf_folderpath, F_OK) == 0)
        {
            cout << subiniconf_folder << " folder exists" << endl;
        }
        else
        {
            folder.createsubIniConf(subiniconf_dir);
        }
        string newconf_file = "ini_conf";
        string newconf_path = subiniconf_dir + "/" + newconf_file;
        string eleconf_file = "ele_conf";
        string eleconf_path = root_dir + "/" + eleconf_file;
        write_poscar(iniconfrlx_path, eleconf_path, newconf_path, n_iter);
        ofstream ofs("status.log",ios::app);
        ofs << n_iter << "  YES  " << endl;
        ofs.clear();
        verlet.conf_md(root_dir,subiniconf_dir,num_step,dump_freq,temp);

    }
    event_num.clear();
}

void KMC::AddEle::write_poscar(string iniconfrlx_path, string eleconf_path, string newconf_path, int n_iter)
{
    opt_cell_ini = io_obj.read_stru(iniconfrlx_path);
    string head = io_obj.header;
    prec_t scale_factor = io_obj.scaling;
    vector<string> elem = io_obj.elements;
    int natom_ini = opt_cell_ini.get_atom_num();
    array<double,9> box = opt_cell_ini.get_box();
    vector<index_t> atype_ini = opt_cell_ini.get_atype();
    vector<prec_t> ini_coords = opt_cell_ini.positions.get_coords();
    opt_cell_ele = io_obj.read_stru(eleconf_path);
    int natom_ele = opt_cell_ele.get_atom_num();
    vector<index_t> atype_ele = opt_cell_ele.get_atype();
    vector<prec_t> ele_coords = opt_cell_ele.positions.get_coords();
    vector<double> z_ini_coords,z_ele_coords;
    z_ini_coords.resize(natom_ini), z_ele_coords.resize(natom_ele);
    for (int i = 0; i < natom_ini; i++)
    {   
        z_ini_coords[i] = ini_coords[3*i+2];
    }   
    for (int i = 0; i < natom_ele; i++)
    {
        z_ele_coords[i] = ele_coords[3*i+2];
    } 

    auto maxElement = max_element(z_ini_coords.begin(), z_ini_coords.end());
    double max_z = *maxElement;
    auto minElement = min_element(z_ele_coords.begin(), z_ele_coords.end());
    double min_z = *minElement;
    double z_diff = max_z - min_z + 1.0;
    ofstream ofs(newconf_path);
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
    unordered_map<index_t, size_t> elenum_ini;
    unordered_map<index_t, size_t> elenum_ele;
    for (int i = 0 ; i < elem.size() ; i++)
        elenum_ini.insert(pair<index_t, size_t>(i, 0)),
        elenum_ele.insert(pair<index_t, size_t>(i, 0));
    for (auto type : atype_ini)
        elenum_ini[type] += 1;
    for (auto type : atype_ele)
        elenum_ele[type] += 1;
    for (int i = 0 ; i < elem.size() ; i++)
        ofs << elenum_ini[i] + elenum_ele[i] << " ";
    ofs << "\n";
    ofs << "Selective Dynamics\n";
    ofs << "Cartesian\n";
    int num_ini = 0, num_ele = 0;
    for (int i = 0 ; i < elem.size() ; i++)
    {
        for (int j = 0 ; j < elenum_ini[i]; j++)
        {
	    num_ini += 1;
	    if (z_ini_coords[num_ini-1] <= (max_z - 3.0))
                ofs << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)] << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)+1] << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)+2] << "    F" << "  F" << "  F" << endl;
            else
                ofs << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)] << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)+1] << setw(22) << setprecision(16) << fixed << ini_coords[3*(num_ini-1)+2] << "    T" << "  T" << "  T" << endl;
        }
        for (int j = 0 ; j < elenum_ele[i]; j++)
        {
	    num_ele += 1;
            ofs << setw(22) << setprecision(16) << fixed << ele_coords[3*(num_ele-1)] << setw(22) << setprecision(16) << fixed << ele_coords[3*(num_ele-1)+1] << setw(22) << setprecision(16) << fixed << ele_coords[3*(num_ele-1)+2]+z_diff << "    T" << "  T" << "  T" << endl;
        }
    }
}
