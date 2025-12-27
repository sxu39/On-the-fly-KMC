#include "conf_gen.h"
#include "fileio.h"
#include "copy.h"
#include "para.h"
#include "folder.h"
#include "opt.h"
#include "neb.h"
#include "model_devi.h"
#include "energy_screen.h"
#include "conf_screen.h"
#include "neb_analysis.h"
#include "fix_atom.h"
#include "add_ele.h"
#include "verlet.h"
#include "collect_traj.h"
#include <dirent.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <filesystem>

using namespace std;
using json = nlohmann::json;

int main(int argc, char *argv[])
{

    ifstream ifs("para.json");
    if (!ifs.is_open()) 
    {
        cerr << "Failed to open para.json" << endl;
        return 1;
    }
    json data;
    try 
    {
        ifs >> data;
    } 
    catch (const json::parse_error& e) 
    {
        cerr << "Failed to parse para.json: " << e.what() << endl;
        return 1;
    }
    string ini_opt = data.value("ini_opt", "YES");
    int iter_num = data.value("iter_num", 1);
    int iter_before = data.value("iter_before", 0);
    int conf_num = data.value("conf_num", 50);
    int cont_num = data.value("cont_num", 15);
    int image_num = data.value("image_num",5);
    int max_iter = data.value("max_iter",100);
    int num_step = data.value("num_step",10000);
    int dump_freq = data.value("dump_freq",1);
    int iter_num_min = data.value("iter_num_min",100);
    int upper_ele_layer = data.value("upper_ele_layer",1);
    int min_event = data.value("min_event",100);
    int max_event = data.value("max_event",500);
    double aveforce_limit = data.value("aveforce_limit",0.1);
    double etol = data.value("etol",1.0e-3);
    double dis_threshold = data.value("dis_threshold",0.5);
    double dis_devi = data.value("dis_devi",0.5);
    double rcut = data.value("rcut",2.0);
    double rcut_fix = data.value("rcut_fix",4.0);
    double radius = data.value("radius",3.0);
    double z_fix = data.value("z_fix",6.0);
    double temp = data.value("temp", 300.0);
    double lower_barrier = data.value("lower_barrier", 0.0);
    double upper_barrier = data.value("upper_barrier", 1.0);
    string dump_traj = data.value("dump_traj", "NO");
    string opt_method = data.value("opt_method","CG");
    string coll_traj = data.value("coll_traj","NO");
    KMC::in_out io_obj;
    KMC::timer timer_obj;
    KMC::neb temp_neb;
    temp_neb.set_timer(&timer_obj);
    temp_neb.set_in_out(&io_obj);
    temp_neb.set_potential("graph.000.pb");
    KMC::opt temp_opt(etol, 0);
    io_obj.set_log("LOG");
    io_obj.set_phy("ENERGY");
    temp_opt.init("graph.000.pb", opt_method, &io_obj, &timer_obj);
    vector<string> modelPaths = {"graph.000.pb","graph.001.pb","graph.002.pb","graph.003.pb"};
    KMC::ModelDevi model(modelPaths);
    Copy copy;
    FileIO fileio;
    Folder folder;
    EnergyScreen energyscreen;
    KMC::FixAtom fixatom;
    NEBAna nebana;
    int natoms;
    string root_dir = getcwd(NULL,0);
    string iniconf_folder = "IniConf";
    string iniconf_dir = root_dir + "/" + iniconf_folder;
    const char* iniconf_folderpath = iniconf_dir.c_str();
    string iniconfrlx_folder = "IniConfRlx";
    string iniconfrlx_dir = root_dir + "/" + iniconfrlx_folder;
    const char* iniconfrlx_folderpath = iniconfrlx_dir.c_str();	
    int num_ele_layer;
    string subiniconf_folder,subiniconf_dir;
    if (ini_opt == "YES")
    {	    
        //create the iniconf folder
        if (access(iniconf_folderpath, F_OK) == 0) 
        {
            cout << iniconf_folder << " folder exists" << endl;
        } 
        else 
        {
            folder.createIniConf(iniconf_dir);
        }
        int num_ele_layer = 1;
        subiniconf_folder = to_string(num_ele_layer) + "Layer";
        subiniconf_dir = iniconf_dir + "/" + subiniconf_folder;
        const char* subiniconf_folderpath = subiniconf_dir.c_str();
        //create the subiniconf folder
        if (access(subiniconf_folderpath, F_OK) == 0)
        {
            cout << subiniconf_folder << " folder exists" << endl;
        }
        else
        {
            folder.createsubIniConf(subiniconf_dir);
        }
        //create the iniconfrlx folder
        if (access(iniconfrlx_folderpath, F_OK) == 0)  
        {    
            cout << iniconfrlx_folder << " folder exists" << endl;
        }    
        else 
        {    
            folder.createIniConfRlx(iniconfrlx_dir);
        }   
        string iniconf_src_path = root_dir + "/ini_conf";
        string iniconf_dst_path = subiniconf_dir + "/ini_conf";
        copy.copy_file(iniconf_src_path, iniconf_dst_path);
        string log_src_path = subiniconf_dir + "/energy.log";
        string log_dst_path = iniconfrlx_dir + "/energy.log";
        string outlog_path = subiniconf_dir + "/out.log";
        string enelog_path = log_src_path;
	io_obj.set_log(outlog_path);
        io_obj.set_phy(enelog_path);
        KMC::cell opt_cell = io_obj.read_stru(iniconf_dst_path);
        temp_opt.set_system(&opt_cell);
        temp_opt.run(1000);
        string iniconfrlx_src_path = subiniconf_dir + "/ini_conf_rlx";
        string iniconfrlx_dst_path = iniconfrlx_dir + "/ini_conf_rlx";
        io_obj.print_stru(opt_cell, iniconfrlx_src_path);
        copy.copy_file(log_src_path, log_dst_path);
        copy.copy_file(iniconfrlx_src_path, iniconfrlx_dst_path);
    }
    ofstream ofs_time("running_time.log", ios::app);
    for (int n_iter = 1 + iter_before; n_iter <= iter_num; n_iter++)
    {
	auto start = std::chrono::high_resolution_clock::now();
        string iniconfrlx_src_path = iniconfrlx_dir + "/ini_conf_rlx";
        string iniconfrlx_dst_path = iniconfrlx_dir + "/ini_conf_rlx_" + to_string(n_iter-1);
        copy.copy_file(iniconfrlx_src_path, iniconfrlx_dst_path);
        string log_src_path = iniconfrlx_dir + "/energy.log";
        string log_dst_path = iniconfrlx_dir + "/energy.log_" + to_string(n_iter-1);
        copy.copy_file(log_src_path, log_dst_path);
        string status_file = "status.log";
        string status_path = root_dir + "/" + status_file;
        Para para_yes = fileio.read_status(status_path);
        int num_yes = para_yes.get_yes_count();
        ostringstream oss;
        oss << "iter_" << setfill('0') << setw(4) << n_iter;
        string iter_name = oss.str();
        string iter_dir = root_dir + "/" + iter_name; 
        const char* iter_folderpath = iter_dir.c_str();
        if (access(iter_folderpath, F_OK) == 0)  
        {   
            cout << iter_name  << " folder exists" << endl;
        }   
        else 
        {   
            folder.createIterFolder(iter_dir);
        }
        string gen_folder = "1.conf_gen";
        string gen_folder_dir = root_dir + "/" + iter_name + "/" + gen_folder; 
        const char* gen_folderpath = gen_folder_dir.c_str();
        string opt_folder = "2.conf_opt";
        string opt_folder_dir = root_dir + "/" + iter_name + "/" + opt_folder;
        const char* opt_folderpath = opt_folder_dir.c_str();
        string confcand_folder = "Candidate";
        string confcand_folder_dir = opt_folder_dir + "/" + confcand_folder; 
        const char* confcand_folderpath = confcand_folder_dir.c_str();
        string conf_folder = "3.conf_candidate";
        string conf_folder_dir = root_dir + "/" + iter_name + "/" + conf_folder;
        const char* conf_folderpath = conf_folder_dir.c_str();
        string neb_folder = "4.conf_neb_cal";
        string neb_folder_dir = root_dir + "/" + iter_name + "/" + neb_folder;
        const char* neb_folderpath = neb_folder_dir.c_str();
        //create the gen_folder and generate the random configurations
        if (access(gen_folderpath, F_OK) == 0)  
        {
            cout << gen_folder << " folder exists" << endl;
            string gen_folder_src = gen_folder_dir;
            string gen_folder_dst = gen_folder_dir + ".bak";
            copy.copy_file(gen_folder_src, gen_folder_dst);  
        }
        else 
        {
            folder.createGenFolder(gen_folder_dir);
            KMC::ConfGen cg;
            cg.genconf(root_dir,gen_folder_dir,iniconfrlx_src_path,rcut,rcut_fix,radius,z_fix,n_iter,num_yes,conf_num);
        }        
        //create the opt_folder and run the optimization
        if (access(opt_folderpath, F_OK) == 0)
        {
            cout << opt_folder << " folder exists" << endl;
            string opt_folder_src = opt_folder_dir;
            string opt_folder_dst = opt_folder_dir + ".bak";
            copy.copy_file(opt_folder_src, opt_folder_dst);                  
        }
        else
        {
            folder.createOptFolder(opt_folder_dir);
	    folder.createConfCandFolder(confcand_folder_dir);
        }	
        //auto opt_start = std::chrono::high_resolution_clock::now();
        for (int n_conf = 1; n_conf <= conf_num; n_conf++)
        {
            string inpos_path = gen_folder_dir + "/trial_conf_" + to_string(n_conf);
            string outpos_path = opt_folder_dir + "/final_conf_rlx_" + to_string(n_conf);	
    	    string outlog_path = opt_folder_dir + "/out.log";
    	    string enelog_path = opt_folder_dir + "/energy.log_" + to_string(n_conf);
	    io_obj.set_log(outlog_path);
            io_obj.set_phy(enelog_path);
            KMC::cell opt_cell = io_obj.read_stru(inpos_path);
	    temp_opt.set_system(&opt_cell);
	    if (dump_traj == "YES")
            {		    
    	        string purt_folder = "purt_" + to_string(n_conf);
                string purt_folder_dir = opt_folder_dir + "/" + purt_folder;
    	        string purt_file = purt_folder_dir +  "/POSCAR_OPT";
                folder.createIterFolder(purt_folder_dir); 
                temp_opt.run(max_iter, dump_freq, purt_file);
	    }
	    else
	    {
                temp_opt.run(max_iter);
	    }
    	    io_obj.print_stru(opt_cell,outpos_path);
        }
        string upperpos_file = "upper_pos.log";
	string upperpos_path = root_dir + "/" + upperpos_file;
        ifstream ifs_pos(upperpos_path, ios::in);
        string pos_line;
        string lastLine;
        while (getline(ifs_pos, pos_line))
        {
            lastLine = pos_line;
        }
        istringstream iss(lastLine);
        int iternum;
	double value1,value2;
        iss >> iternum >> value1 >> value2;
        double z_low = value1;
        double z_high = value2;
	double upper_devi;
	string modeldevi_file = "upper_devi.log";
        string modeldevi_path = root_dir + "/" + modeldevi_file;
        ofstream ofs_devi(modeldevi_path,ios::app); 
        for (int i = n_iter - 1; i <= n_iter - 1; i++)
	{       	
            upper_devi = model.cal_devi_s(iniconfrlx_dst_path,z_low,z_high);
            ofs_devi << n_iter - 1 << " " << upper_devi << endl;
	}     
        model.cal_devi_m(opt_folder_dir,conf_num,z_low,z_high);
        energyscreen.energy_screen(root_dir,iniconfrlx_dir,opt_folder_dir,confcand_folder_dir,conf_num,max_iter,upper_devi,aveforce_limit);

        //create the screened conf folder
        if (access(conf_folderpath, F_OK) == 0)
        {
            cout << conf_folder << " folder exists" << endl;
            string conf_folder_src = conf_folder_dir;
            string conf_folder_dst = conf_folder_dir + ".bak";
            copy.copy_file(conf_folder_src, conf_folder_dst);
        }
        else
        {
            folder.createConfFolder(conf_folder_dir);
        }
	KMC::ConfScreen cs;
        cs.conf_screen(iniconfrlx_dir, confcand_folder_dir, conf_folder_dir, dis_threshold);

        //create the neb folder
        if (access(neb_folderpath, F_OK) == 0)
        {
            cout << neb_folder << " folder exists" << endl;
            string neb_folder_src = neb_folder_dir;
            string neb_folder_dst = neb_folder_dir + ".bak";
            copy.copy_file(neb_folder_src, neb_folder_dst);
        }
        else
        {
            folder.createNEBFolder(neb_folder_dir);
        }
     	
        string src_folder_path = conf_folder_dir;
        string dst_folder_path = neb_folder_dir;
        string cmd_copy = "cp -rf " + src_folder_path + "/* " + dst_folder_path;
        system(cmd_copy.c_str());
        string folder_path = dst_folder_path;
        int n_folder = 0;
        //auto neb_start = std::chrono::high_resolution_clock::now();
        if (chdir(folder_path.c_str()) == 0)
        {
            DIR* dir = opendir(".");
            if (dir != NULL)
            {
                struct dirent* entry;
                while ((entry = readdir(dir)) != NULL)
                {
                    string file_name = entry->d_name;
                    if (file_name.substr(0, 14) == "final_conf_rlx")
                    {
                        string newFolderName = "1_" + file_name;
                        if (mkdir(newFolderName.c_str(), 0777) == 0)
                        {
                            n_folder += 1;
                            rename(file_name.c_str(), (newFolderName + "/final_conf_rlx").c_str());
                            string subfolder_path =  newFolderName;
                            string finalconf_file = "final_conf_rlx";
                         
                            // read initial and final structures, and add them to neb module
               		    string log_path = neb_folder_dir + "/" + subfolder_path + "/out.log";
                            string energy_path = neb_folder_dir + "/" + subfolder_path + "/energy.log";
                            io_obj.set_log(log_path);
                            io_obj.set_phy(energy_path);

            	  	    string initial_file_path = iniconfrlx_src_path;
            		    string final_file_path = neb_folder_dir + "/" + subfolder_path + "/" + finalconf_file;
                            //string newinitial_file_path = neb_folder_dir + "/" + subfolder_path + "/ini_conf_rlx_new";
                            //string newfinal_file_path = neb_folder_dir + "/" + subfolder_path + "/final_conf_rlx_new";
			    string newinitial_file_path = neb_folder_dir + "/" + subfolder_path + "/POSCAR_0";
                            string newfinal_file_path = neb_folder_dir + "/" + subfolder_path + "/POSCAR_"+to_string(image_num-1);
            		    fixatom.fix_atom(initial_file_path,final_file_path,newinitial_file_path,newfinal_file_path,dis_devi);
                            KMC::cell initial_cell = io_obj.read_stru(newinitial_file_path);
                            KMC::cell final_cell = io_obj.read_stru(newfinal_file_path);
                            temp_neb.set_system(&initial_cell, &final_cell);
                            temp_neb.interpolate();
            		    fixatom.fix_image(initial_file_path,final_file_path,newinitial_file_path,newfinal_file_path);
            		    KMC::cell initial_cell_new = io_obj.read_stru(newinitial_file_path);
                            KMC::cell final_cell_new = io_obj.read_stru(newfinal_file_path);
                            temp_neb.set_system(&initial_cell_new, &final_cell_new);
                            // perform neb minimization
                            temp_neb.run();
            		    string temp_pos_path =  neb_folder_dir + "/" + subfolder_path;
                            temp_neb.out_stru(temp_pos_path);
                            //temp_neb.out_stru();
            		    string nebfile_name = "neb_ene.txt";
            		    string nebfile_src_path = neb_folder_dir + "/" + nebfile_name;
            		    string nebfile_dst_path = neb_folder_dir + "/" + subfolder_path + "/" + nebfile_name;
            		    copy.move_file(nebfile_src_path, nebfile_dst_path);
            	        }
                    }
                }
                //closedir(dir);
            }
            closedir(dir);
        }
	nebana.neb_analysis(root_dir,iniconfrlx_dir,opt_folder_dir,neb_folder_dir,n_iter,image_num,temp,lower_barrier,upper_barrier);
	int iter_num_limit = 0, last_integer;
        ifstream ifs(status_path);
        if (ifs.peek() == ifstream::traits_type::eof()) 
        {
            iter_num_limit = iter_num_min; 
        }
        else 
        {
            string buff;
            vector<string> lines;
            while (getline(ifs, buff)) 
            {
                lines.push_back(buff);
            }
            //ifs.close();
            if (!lines.empty()) 
            {
                string last_line = lines.back();
                istringstream ss(last_line);
                string str;
                ss >> last_integer >> str;
                iter_num_limit = last_integer + iter_num_min;
		//cout << "Last Integer is   " << last_integer << endl;
            }
        }
	ifs.close();
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = end - start;
        ofs_time << n_iter << " " << duration.count() << endl;
	//count the effective KMC number
        string iter_file = "iter.log";
        string iter_path = root_dir + "/" + iter_file;
        ifstream iter_log(iter_path,ios::in);
        string line,value;
	int eventnum;
	vector<int> eff_num(10, 0);
        if (num_yes == 0)
	{	
            for (int iter = 1; iter <= n_iter; iter++)
            {
                for (int j = 0; j < 8; j++)
                {
                    getline(iter_log,line);
                    if (line.find("Event_num") != string::npos)
                    {
                        stringstream ss(line);
                        ss >> value >> eventnum;
                        if (eventnum > 0)
                        {	
                            eff_num[num_yes] += 1; 
			    //cout << n_iter << " " << num_yes << " " << eff_num[num_yes] << endl;
                        }	    
                    }
                }
	    }
	    cout << n_iter << " " << num_yes << " " << eff_num[num_yes] << endl;
	}
	else if(num_yes > 0)
	{
	    for (int iter = 1; iter <= last_integer; iter++)
            {
                for (int j = 0; j < 8; j++)
                {
                    getline(iter_log,line);
		}
            }
	    for (int iter = last_integer + 1; iter <= n_iter; iter++)
            {
                for (int j = 0; j < 8; j++)
                {
                    getline(iter_log,line);
                    if (line.find("Event_num") != string::npos)
                    {
                        stringstream ss(line);
                        ss >> value >> eventnum;
                        if (eventnum > 0)
                        {
                            eff_num[num_yes] += 1;
			    //cout << n_iter << " " << num_yes << " " << eff_num[num_yes] << endl;
                        }
                    }
                }
            }
	    cout << n_iter << " " << num_yes << " " << eff_num[num_yes] << endl;
        }
	//add the electrolyte and do the NVE simulation
	if (n_iter > iter_num_limit && eff_num[num_yes] >= min_event)
	{	
	    num_ele_layer = num_yes + 2;
            KMC::AddEle ae;
            ae.addele(root_dir,iniconf_dir,iniconfrlx_dir,n_iter,cont_num,num_yes,num_step,dump_freq,max_event,temp);
            subiniconf_folder = to_string(num_ele_layer)+ "Layer";
            subiniconf_dir = iniconf_dir + "/" + subiniconf_folder;
            auto folderExists = [](const string& subiniconf_dir) 
            {
                struct stat info;
                if (stat(subiniconf_dir.c_str(), &info) != 0) 
                {
                    return false;
                }
                return (info.st_mode & S_IFDIR) != 0;
            };
            if (folderExists(subiniconf_dir)) 
            {
                string conf_file = "POSCAR_" + to_string(num_step);
                string conf_path = subiniconf_dir + "/" + conf_file;
                string newconf_file = "ini_conf_new";
                string newconf_path = subiniconf_dir + "/" + newconf_file;
                fixatom.fix_height(conf_path,newconf_path,z_fix);
                string log_src_path = subiniconf_dir + "/energy.log";
                string log_dst_path = iniconfrlx_dir + "/energy.log";
                string outlog_path = subiniconf_dir + "/out.log";
                string enelog_path = log_src_path;
                io_obj.set_log(outlog_path);
                io_obj.set_phy(enelog_path);
                KMC::cell opt_cell = io_obj.read_stru(newconf_path);
                temp_opt.set_system(&opt_cell);
                temp_opt.run(1000);
                iniconfrlx_src_path = subiniconf_dir + "/ini_conf_rlx";
                iniconfrlx_dst_path = iniconfrlx_dir + "/ini_conf_rlx";
                io_obj.print_stru(opt_cell, iniconfrlx_src_path);
                copy.copy_file(log_src_path, log_dst_path);
                copy.copy_file(iniconfrlx_src_path, iniconfrlx_dst_path);
            }
            Para para = fileio.read_status(status_path);
            int num_yes_new = para.get_yes_count();
	    if (num_yes_new > upper_ele_layer - 1)
                break;
	}
    }
    if (coll_traj == "YES" or coll_traj == "yes" or coll_traj == "Yes")
    {    
        KMC::CollTraj cj;
        cj.colltraj(root_dir,iter_num,conf_num);
    }    
    return 0;   
}

