#include "model_devi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

KMC::ModelDevi::ModelDevi(vector<string>& modelPaths)
{
    for (const std::string& path : modelPaths)
    {
        deepmd::DeepPot* dp = new deepmd::DeepPot();
        dp->init(path);
        dp_models.push_back(dp);
    }
}

void KMC::ModelDevi::cal_devi_m(string opt_folder_dir,int conf_num, double z_low, double z_high)
{
   for (int n_conf = 1; n_conf <= conf_num; n_conf++)
   {
       string conf_file = "final_conf_rlx_" + to_string(n_conf);
       string conf_path = opt_folder_dir + "/" + conf_file; 
       opt_cell = io_obj.read_stru(conf_path);
       int natoms = opt_cell.get_atom_num();
       array<double,9> box = opt_cell.get_box();
       vector<index_t> atype = opt_cell.get_atype();
       vector<prec_t> coords = opt_cell.positions.get_coords();
       double e;
       vector <double> f,v,f_ave_x,f_ave_y,f_ave_z,f_devi;
       vector <vector<double>> fx(4,vector<double>()),fy(4,vector<double>()),fz(4,vector<double>());
       f_ave_x.resize(natoms),f_ave_y.resize(natoms),f_ave_z.resize(natoms),f_devi.resize(natoms);
       vector<double> boxs;
       boxs.clear();
       f_devi.clear();
       for (int i = 0 ; i < 9 ; i++)
           boxs.push_back(box[i]);
       for (int i = 0; i < dp_models.size(); i++)
       {
           deepmd::DeepPot* dp = dp_models[i];
           dp->compute(e, f, v, coords, atype, boxs);
           for (int n = 0; n < natoms; n++)
           {
               fx[i].push_back(f[3*n]);
               fy[i].push_back(f[3*n+1]);
               fz[i].push_back(f[3*n+2]);
           }
       } 
       for (int n = 0; n < natoms; n++)
       {
           f_ave_x.clear();
           f_ave_y.clear();
           f_ave_z.clear();
	   if (coords[3*n+2] > z_low && coords[3*n+2] < z_high)
           {
               for (int i = 0; i < 4; i++)
               {
                   f_ave_x[n] += fx[i][n];
                   f_ave_y[n] += fy[i][n];
                   f_ave_z[n] += fz[i][n];
               }    
               f_ave_x[n] = f_ave_x[n]/4;
               f_ave_y[n] = f_ave_y[n]/4;
               f_ave_z[n] = f_ave_z[n]/4;
               double aaa = 0.0;    
               for (int i = 0; i < 4; i++)
               {
                   double devi = (fx[i][n] - f_ave_x[n])*(fx[i][n] - f_ave_x[n]) + (fy[i][n] - f_ave_y[n])*(fy[i][n] - f_ave_y[n]) + (fz[i][n] - f_ave_z[n])*(fz[i][n] - f_ave_z[n]);
                   aaa += devi;
               }
               double bbb = sqrt(aaa/4);
               f_devi.push_back(bbb);
       	   }
       }
       double max_force_devi = *std::max_element(f_devi.begin(), f_devi.end());
       double min_force_devi = *std::min_element(f_devi.begin(), f_devi.end());
       double ave_force_devi = 0.0;
       for (int i = 0; i < f_devi.size(); i++)
       {
           ave_force_devi += f_devi[i];
       }
       ave_force_devi = ave_force_devi / f_devi.size();
       string devi_log_file = "model_devi.log";
       string devi_log_path = opt_folder_dir + "/" + devi_log_file;
       ofstream ofs(devi_log_path,ios::app);
       ofs << n_conf << " " << max_force_devi << " " << min_force_devi << " " << ave_force_devi << endl;
   }
}

double KMC::ModelDevi::cal_devi_s(string filename, double z_low, double z_high)
{
   opt_cell = io_obj.read_stru(filename);
   int natoms = opt_cell.get_atom_num();
   array<double,9> box = opt_cell.get_box();
   vector<index_t> atype = opt_cell.get_atype();
   vector<prec_t> coords = opt_cell.positions.get_coords();
   double e;
   vector <double> f,v,f_ave_x,f_ave_y,f_ave_z,f_devi;
   vector <vector<double>> fx(4,vector<double>()),fy(4,vector<double>()),fz(4,vector<double>());
   f_ave_x.resize(natoms),f_ave_y.resize(natoms),f_ave_z.resize(natoms),f_devi.resize(natoms);
   vector<double> boxs;
   boxs.clear();
   f_devi.clear();
   for (int i = 0 ; i < 9 ; i++)
       boxs.push_back(box[i]);
   for (int i = 0; i < dp_models.size(); i++)
   {
       deepmd::DeepPot* dp = dp_models[i];
       dp->compute(e, f, v, coords, atype, boxs);
       for (int n = 0; n < natoms; n++)
       {
           fx[i].push_back(f[3*n]);
           fy[i].push_back(f[3*n+1]);
           fz[i].push_back(f[3*n+2]);
       }
   }
   for (int n = 0; n < natoms; n++)
   {
       f_ave_x.clear();
       f_ave_y.clear();
       f_ave_z.clear();
       if (coords[3*n+2] > z_low && coords[3*n+2] < z_high)
       {
           for (int i = 0; i < 4; i++)
           {
               f_ave_x[n] += fx[i][n];
               f_ave_y[n] += fy[i][n];
               f_ave_z[n] += fz[i][n];
           }
           f_ave_x[n] = f_ave_x[n]/4;
           f_ave_y[n] = f_ave_y[n]/4;
           f_ave_z[n] = f_ave_z[n]/4;
           double aaa = 0.0;
           for (int i = 0; i < 4; i++)
           {
               double devi = (fx[i][n] - f_ave_x[n])*(fx[i][n] - f_ave_x[n]) + (fy[i][n] - f_ave_y[n])*(fy[i][n] - f_ave_y[n]) + (fz[i][n] - f_ave_z[n])*(fz[i][n] - f_ave_z[n]);
               aaa += devi;
           }
           double bbb = sqrt(aaa/4);
           f_devi.push_back(bbb);
       }
   }
   double max_force_devi = *std::max_element(f_devi.begin(), f_devi.end());
   double min_force_devi = *std::min_element(f_devi.begin(), f_devi.end());
   double ave_force_devi = 0.0;
   for (int i = 0; i < f_devi.size(); i++)
   {
       ave_force_devi += f_devi[i];
   }
   ave_force_devi = ave_force_devi / f_devi.size();
   return ave_force_devi;
}
