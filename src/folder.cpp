#include "folder.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>

using namespace std;

string confcand_folder_dir;
string conf_folder_dir;
string neb_folder_dir;

void Folder::createEleConf(string eleconf_dir)
{
    string cmd_mkdir = "mkdir " + eleconf_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createIniConf(string iniconf_dir)
{
    string cmd_mkdir = "mkdir " + iniconf_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createsubIniConf(string subiniconf_dir)
{
    string cmd_mkdir = "mkdir " + subiniconf_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createIniConfRlx(string iniconfrlx_dir)
{
    string cmd_mkdir = "mkdir " + iniconfrlx_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createIterFolder(string iter_dir)
{
    string cmd_mkdir = "mkdir " + iter_dir;
    (void)system(cmd_mkdir.c_str());
}

void Folder::createGenFolder(string gen_folder_dir)
{
    string cmd_mkdir = "mkdir " + gen_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createOptFolder(string opt_folder_dir)
{
    string cmd_mkdir = "mkdir " + opt_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createConfCandFolder(string confcand_folder_dir)
{
    string cmd_mkdir = "mkdir " + confcand_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createConfFolder(string conf_folder_dir)
{
    string cmd_mkdir = "mkdir " + conf_folder_dir;
    system(cmd_mkdir.c_str()); 
}

void Folder::createNEBFolder(string neb_folder_dir)
{
    string cmd_mkdir = "mkdir " + neb_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createPurtFolder(string purt_folder_dir)
{
    string cmd_mkdir = "mkdir " + purt_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createTaskFolder(string task_folder_dir)
{
    string cmd_mkdir = "mkdir " + task_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createModeldeviFolder(string modeldevi_folder_dir)
{
    string cmd_mkdir = "mkdir " + modeldevi_folder_dir;
    system(cmd_mkdir.c_str());
}

void Folder::createTrajFolder(string traj_folder_dir)
{
    string cmd_mkdir = "mkdir " + traj_folder_dir;
    system(cmd_mkdir.c_str());
}
