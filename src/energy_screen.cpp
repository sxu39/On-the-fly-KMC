#include "energy_screen.h"
#include "folder.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <algorithm>

using namespace std;

void EnergyScreen::energy_screen(string root_dir, string iniconfrlx_dir, string opt_folder_dir, string confcand_folder_dir, int conf_num, int max_iter, double upper_devi, double aveforce_limit)
{
   int step;
   double E0, energy;
   double max_devi_v, min_devi_v, avg_devi_v, max_devi_f, min_devi_f, avg_devi_f;
   string devi_file = "model_devi.log";
   string devi_path = opt_folder_dir + "/" + devi_file;
   double max_force_devi[conf_num],min_force_devi[conf_num],ave_force_devi[conf_num];
   int num;
   ifstream ifs(devi_path,ios::in);
   for (int i = 1; i <= conf_num; i++)
   {
       ifs >> num >> max_force_devi[i] >> min_force_devi[i] >> ave_force_devi[i];
   } 
   string inilog_file = "energy.log";
   string inilog_path = iniconfrlx_dir + "/" + inilog_file;
   extract_optenergy(inilog_path, energy);
   E0 = energy;
   string energy_file = "ene_dp.txt";
   string energy_path = opt_folder_dir + "/" + energy_file;
   ofstream ofs(energy_path,ios::out);
   ofs << "0" << " " << fixed << setprecision(4) << energy << endl;
   upper_devi = min(upper_devi + 0.01, aveforce_limit);
   int n_cand = 0;
   string cand_file = "num_cand.log";
   string cand_path = root_dir + "/" + cand_file;
   ofstream ofs_cand(cand_path,ios::app);
   for (int i = 1; i <= conf_num; i++)
   {
       string log_file = "energy.log_" + to_string(i);
       string log_path = opt_folder_dir + "/" + log_file;
       int line_num = 0;
       string line;
       ifstream ifs(log_path);
       while (getline(ifs, line))
       {
           line_num++;
       }
       extract_optenergy(log_path, energy);
       double energy_limit = E0 + 1.0;
       if (ave_force_devi[i] <= upper_devi)
       //if (energy < energy_limit && (line_num - 3) < max_iter)
       {   
	   n_cand += 1;
	   if (energy < energy_limit && (line_num - 3) < max_iter)
	   //if (ave_force_devi[i] < upper_devi)
           {
               stringstream ss("");
               ss << "final_conf_rlx_" << i;
               string src_path = opt_folder_dir + "/" + ss.str();
               string dst_path = confcand_folder_dir; + "/" + ss.str();
               string cmd_copy = "cp " + src_path + " " + dst_path;
               system(cmd_copy.c_str());
           }
       }           
       ofs << i << " " << fixed << setprecision(4) << energy << endl;
   }
   ofs.close();
   double ratio = static_cast<double>(conf_num - n_cand) / static_cast<double>(conf_num);
   ofs_cand << conf_num - n_cand << " " << ratio << endl; 
}

void EnergyScreen::extract_optenergy(string file_name, double& energy)
{
    ifstream ifs(file_name);
    string line;
    string match = "Final Energy";
    while (getline(ifs, line))
    {
        if (line.find(match) != string::npos)
        {
            size_t pos = line.find_last_of(" ");
            string value = line.substr(pos + 1);
            energy = stod(value);
            return;
        }
    }
    ifs.close();
}

