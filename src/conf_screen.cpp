#include "conf_screen.h"
#include "folder.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <vector>
#include <unistd.h>
#include <cstdlib>
#include <dirent.h>
#include <cstring>

using namespace std;

void KMC::ConfScreen::conf_screen(string iniconfrlx_dir, string confcand_folder_dir, string conf_folder_dir, double dist_threshold)
{
    string conf_file = "ini_conf_rlx";
    string conf_path = iniconfrlx_dir + "/" + conf_file; 
    opt_cell = io_obj.read_stru(conf_path);
    int natoms = opt_cell.get_atom_num();
    int atomid,atomtype;
    vector<vector<int>> n(natoms, vector<int>(natoms));
    vector<vector<double>> x(natoms, vector<double>(natoms));
    vector<vector<double>> y(natoms, vector<double>(natoms));
    vector<vector<double>> z(natoms, vector<double>(natoms));
    double x1, y1, z1;
    struct dirent *ptr;
    DIR *dir;
    string line;
    vector <string> files;
    dir = opendir(confcand_folder_dir.c_str());
    files.clear();
    while ((ptr = readdir(dir)) != NULL)
    {
        if (ptr->d_name[0] == '.')  
        {
            continue;
        }
        if (strncmp(ptr->d_name, "final", 5) == 0) 
        {
            files.push_back(ptr->d_name);
        }
    }   
    int aa = files.size();
    for (int i = 0; i < files.size(); i++)
    {  
	conf_path = confcand_folder_dir + "/" + files[i];
        ifstream ifs;
        ifs.open(conf_path, ios::in);
        for (int j = 0; j < 9; j++)
        {
            getline(ifs, line);
        }
        for (int j = 0; j < natoms; j++)
        {   
	    string value1,value2,value3;
            ifs >> x1 >> y1 >> z1 >> value1 >> value2 >> value3;
            x[i][j] = x1;
            y[i][j] = y1;
            z[i][j] = z1;
        }   
        ifs.close(); 
    }
    for (int i = 0; i < natoms; i++)
    {
        for (int j = 0; j < natoms; j++)
        {
            n[i][j] = 0;
        } 
    }
    for (int i = 0; i < files.size(); i++)
    {
        for (int j = 0; j < files.size(); j++)
        {
            for (int k = 0; k < natoms; k++)
            {
                double dx = x[i][k] - x[j][k];
                double dy = y[i][k] - y[j][k];
                double dz = z[i][k] - z[j][k];
                if (abs(dx) < dist_threshold && abs(dy) < dist_threshold && abs(dz) < dist_threshold)
                {
                    n[i][j] = n[i][j] + 1;
                }   
            }
        }
    }
    for (int i = 0; i < files.size(); i++)
    {
        for (int j = 0; j < files.size(); j++)
        {
            if (i != j && n[i][j] == natoms)
            {
                files[i] = files[j];
            } 
        }
        string file_name1 = confcand_folder_dir + "/" + files[i];
        string file_name2 = conf_folder_dir + "/" + files[i];
        fstream ifs(file_name1, ios::in);
        fstream ofs(file_name2, ios::out);
        ifs.seekg(0, ios::end);
        int total = ifs.tellg();
        ifs.seekg(0,ios::beg);
        char buf[100];
        int readCount = 0;
        int count = 0;
        while (count < total)
        {
            ifs.read(buf, sizeof(buf));
            readCount = ifs.gcount();
            count += readCount;
            ofs.write(buf, readCount);
        }   
        ifs.close();
        ofs.close();
    }
    files.clear(); 
    closedir(dir);
}

