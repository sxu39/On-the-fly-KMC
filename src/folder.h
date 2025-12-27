#ifndef FOLDER_H
#define FOLDER_H

#include <string>

using namespace std;

class Folder 
{
    public:
        void createEleConf(string eleconf_dir);
        void createIniConf(string iniconf_dir);
        void createsubIniConf(string subiniconf_dir);
        void createIniConfRlx(string iniconfrlx_dir);
        void createIterFolder(string iter_dir);
        void createGenFolder(string gen_folder_dir);
        void createOptFolder(string opt_folder_dir);
        void createConfCandFolder(string confcand_folder_dir);
        void createConfFolder(string conf_folder_dir);
        void createNEBFolder(string neb_folder_dir);
        void createPurtFolder(string purt_folder_dir);
        void createTaskFolder(string task_folder_dir);
        void createModeldeviFolder(string modeldevi_folder_dir);
	void createTrajFolder(string traj_folder_dir);
};
#endif
