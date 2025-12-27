#include "collect_traj.h"
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sstream>
#include <iomanip>

void KMC::CollTraj::colltraj(string root_dir, int iter_num, int conf_num)
{
     string dst_folder = "01.model_devi";
     string dst_folder_dir = root_dir + "/" + dst_folder;
     folder.createModeldeviFolder(dst_folder_dir);

     for (int n_iter = 1; n_iter <= iter_num; n_iter++)
     {
	 ostringstream oss;
         oss << "iter_" << setfill('0') << setw(4) << n_iter;
    	 string iter_name = oss.str();
 	 string opt_folder = "2.conf_opt";
 	 string opt_folder_dir = root_dir + "/" + iter_name + "/" + opt_folder;
         ostringstream os;
         os << "task." << to_string(iter_num) << "." << setfill('0') << setw(6) << n_iter-1;
         string task_folder = os.str();
         string task_folder_dir = dst_folder_dir + "/" + task_folder;
         folder.createTaskFolder(task_folder_dir);
         string traj_folder_dir = task_folder_dir + "/traj";
	 folder.createTrajFolder(traj_folder_dir);
         int num_file = 0;
	 for (int n_purt = 1; n_purt <= conf_num; n_purt++)
         {
             string purt_folder = "purt_" + to_string(n_purt);
             string purt_folder_dir = opt_folder_dir + "/" + purt_folder;
             DIR* purt_folder_path = opendir(purt_folder_dir.c_str());
             dirent* fileEntry;
             while ((fileEntry = readdir(purt_folder_path)) != nullptr)
             {
                string fileName = fileEntry->d_name;
                //if (fileName.substr(fileName.find_last_of('.') + 1) == "lammpstrj")
		if (fileName.substr(0, 10) == "POSCAR_OPT")
                {
                    string src_file = purt_folder_dir + "/" + fileName;
                    string dst_file = traj_folder_dir + "/" + to_string(num_file) + ".poscar";
		    ifstream ifs_src(src_file, ios::binary);
		    ofstream ofs_dst(dst_file, ios::binary);
                    ofs_dst << ifs_src.rdbuf();
		    string dump_file = traj_folder_dir + "/" + to_string(num_file) + ".lammpstrj";
                    format.pos2dump(dst_file,dump_file);
                    num_file++;
                }
            }
            closedir(purt_folder_path);
        }
    }
}

