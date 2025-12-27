#include <iostream>
#include <fstream>
#include <iomanip> 
#include <random>
#include <unistd.h>
#include <dirent.h>
#include <algorithm>
#include "neb_analysis.h"

const double planck_const = 6.626e-34; //unit: J.S
const double boltz_const = 1.38e-23; //unit: J/K
const double unit_factor = 23.0609*4.184*1000; //ev to J/mol
const double avogadro_const = 6.02e23; //unit: mol-1

void NEBAna::neb_analysis(string root_dir,string iniconfrlx_dir,string opt_folder_dir,string neb_folder_dir,int n_iter,int image_num,double temp,double lower_barrier,double upper_barrier)
{
    double rate_sum,reaction_rate_const;
    vector<double> neb_barrier_list,rate_list,neb_energy,barrier_list;
    vector<string> neb_folder_list,folders,ts_file_list,is_file_list;
    string iterlog_file = "iter.log";
    string iterlog_path = root_dir + "/" + iterlog_file;
    ofstream ofs_iter(iterlog_path,ios::app);
    string event_list_file = "event.log";
    string event_list_path = root_dir + "/" + event_list_file;
    ofstream ofs_list(event_list_path,ios::app);
    ofs_iter << left;
    ofs_iter << setw(12) << "# Step_num  " << setw(10) << n_iter << endl;
    ofs_list << left;
    ofs_list << setw(12) << "# Step_num  " << setw(10) << n_iter << endl;
    rate_sum = 0.0;
    string ts_folder = neb_folder_dir + "/transition_state";
    string cmd_mkdir = "mkdir " + ts_folder;
    system(cmd_mkdir.c_str());
    DIR *dir;
    dir = opendir(neb_folder_dir.c_str());
    if (dir == nullptr) 
    {
        cerr << "Failed to open directory." << endl;
    }
    struct dirent *ent;
    while ((ent = readdir(dir)) != nullptr) 
    { 
        string folder_name = ent->d_name; 
        if (folder_name.find("1_final_conf_rlx") != string::npos) 
        { 
            folders.push_back(folder_name);
        }
    }
    closedir(dir);
    int event_num = 0, num_ts = 0;
    for (auto folder : folders) 
    {   
	string nebene_file = "neb_ene.txt";
        string nebene_path = neb_folder_dir + "/" + folder + "/" + nebene_file;
        vector<double> neb_ene_list;
        ifstream ifs(nebene_path,ios::in);
        int image_id;
        double ene[image_num],energy;
        for (int i = 0;  i < image_num; i++)
        {
            ifs >> image_id >> energy;
            neb_ene_list.push_back(energy);
            ifs.ignore(256,'\n');
            ene[i] = energy;
        }
        ifs.close();
        int peak_num = 0;
        //double lower_barrier = 0.1, upper_barrier = 0.9;
        for (int i = 1;  i < image_num - 1; i++)
        {
            double ene_before, ene_curr, ene_next;
            ene_before = ene[i-1], ene_curr = ene[i], ene_next = ene[i+1];
            if (ene_before < ene_curr && ene_next < ene_curr)
            {
                peak_num++;
            }
        }
        if (peak_num >= 1 ) 
	//if (peak_num == 1 )
        {
	    auto max_elem = max_element(neb_ene_list.begin()+1, neb_ene_list.end()-1);
            int ts_index = distance(neb_ene_list.begin(),max_elem);
	    if (neb_ene_list[ts_index] > neb_ene_list[4])
	    {
                if (neb_ene_list[ts_index] > neb_ene_list[ts_index - 1] && neb_ene_list[ts_index] > neb_ene_list[ts_index + 1]) 
	        {    
	            double max_value = *max_elem;
	            auto min_elem = min_element(neb_ene_list.begin(),max_elem);
	            double min_value = *min_elem;
	            int is_index = distance(neb_ene_list.begin(),min_elem);
                    double barrier = max_value - min_value; 
	            barrier_list.push_back(barrier);
                    if ( barrier > lower_barrier && barrier < upper_barrier)
                    {  
	                 num_ts++;
	                 string neb_file = "POSCAR_" + to_string(ts_index);
	                 string low_ene_file = "POSCAR_" + to_string(is_index);
                         string ts_src = neb_folder_dir + "/" + folder + "/" + neb_file;
                         string ts_dst = neb_folder_dir + "/" + folder + "/POSCAR_TS";
	                 string ts_folder_path = ts_folder + "/POSCAR_" + to_string(num_ts);
                         copy.copy_file(ts_src, ts_dst); 
	                 copy.copy_file(ts_src, ts_folder_path);
	                 ts_file_list.push_back(neb_file);
                         is_file_list.push_back(low_ene_file);
	                 event_num++;
                         neb_barrier_list.push_back(barrier);
                         reaction_rate_const = boltz_const*temp / planck_const * exp(-barrier*unit_factor/(avogadro_const*temp*boltz_const));
	                 //cout << "rate constant is " <<reaction_rate_const <<endl;
                         rate_list.push_back(reaction_rate_const); 
                         neb_folder_list.push_back(folder);
                         rate_sum += reaction_rate_const;
	                 //cout << "rate sum is " << rate_sum << endl;
                    }
	        }
	    }
        }
    }
    ofs_list << setw(12) << "Event_num  " << setw(10) << event_num << endl;
    for (int i = 0; i < barrier_list.size(); i++)
    {
        ofs_list << "  Barrier  " << barrier_list[i] << endl;
    }
    for (int i = 0; i < neb_barrier_list.size(); i++)
    {
        ofs_list << "  " << neb_folder_list[i] << "   " << is_file_list[i] << "   " << ts_file_list[i] << "   " << neb_barrier_list[i] << endl;
    }
    ofs_list.close();
    ofs_iter << setw(12) << "Event_num  " << setw(10) << event_num << endl;
    if (event_num == 0)
    {
        ofs_iter << setw(12) << "Rand_num1  " << setw(10) << "NULL" << endl;
        ofs_iter << setw(12) << "Barrier  " << setw(10) << "NULL" << endl;
	ofs_iter << setw(12) << "TS_conf  " << setw(10) << "NULL" << endl;
        ofs_iter << setw(12) << "Rand_num2  " << setw(10) << "NULL" << endl;
        ofs_iter << setw(12) << "Evol_time  " << setw(10) << "NULL" << endl;
        ofs_iter << setw(12) << "Event  " << setw(10) << "NULL" << endl;
    }
    else if (event_num > 0)
    {
        double rand_num1 = randseed.get_rand_seed();
        double rate_accum = 0.0;
        double rate_tmp = rand_num1 * rate_sum;
        ofs_iter << setw(12) << "Rand_num1  " << setw(10) << rand_num1 << endl;
        string new_ini_conf_folder; 
        for(int i = 0; i < rate_list.size(); i++)
        {  
            rate_accum += rate_list[i];
            if(rate_tmp < rate_accum)
            {
                new_ini_conf_folder = neb_folder_list[i];
                ofs_iter << setw(12) << "Barrier  " << setw(10) << neb_barrier_list[i] << endl;
		ofs_iter << setw(12) << "TS_conf  " << setw(10) << ts_file_list[i] << endl;
		//string ts_src = neb_folder_dir + "/" + new_ini_conf_folder + "/" + ts_file_list[i];
		//string ts_dst = neb_folder_dir + "/" + new_ini_conf_folder + "/POSCAR_TS";
		//copy.copy_file(ts_src, ts_dst);
                break;
            }
        }
        ofs_iter.clear();
        double rand_num2 = randseed.get_rand_seed();
        ofs_iter << setw(12) << "Rand_num2  " << setw(10) << rand_num2 << endl;
        double time_step, simulation_time;
        time_step = -log(rand_num2) / rate_sum;
	//cout << "rate sum is " << rate_sum <<endl;
        simulation_time = time_step*1.0e12;
	//cout << "Time is " <<simulation_time <<endl;
        ofs_iter << setw(12) << "Evol_time  " << setw(10) << simulation_time << endl;
        string final_conf_folder = new_ini_conf_folder;
        int pos = final_conf_folder.find_last_not_of("0123456789");
        string str_part = final_conf_folder.substr(0, pos + 1);
        string num_part = final_conf_folder.substr(pos + 1);
        int num = stoi(num_part);
        string inilog_file = "energy.log_" + to_string(num);
        string inilog_path = opt_folder_dir + "/" + inilog_file;
        string event_file = "final_conf_rlx_" + to_string(num);
        ofs_iter << setw(12) << "Event  " << setw(10) << event_file << endl;
        string newlog_file = "energy.log";
        string newlog_path = iniconfrlx_dir + "/" + newlog_file;
        copy.copy_file(inilog_path, newlog_path);
        string final_conf_name = "final_conf_rlx";
        string final_conf_path = neb_folder_dir + "/" + final_conf_folder + "/" + final_conf_name;
        string initial_conf_name = "ini_conf_rlx";
        string initial_conf_path = iniconfrlx_dir + "/" + initial_conf_name;
	copy.copy_file(final_conf_path, initial_conf_path);
    }
    if (chdir("../../") != 0) 
    {
        perror("Failed to change directory");
    }    
}

