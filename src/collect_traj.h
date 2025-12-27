#ifndef COLLECT_TRAJ_H
#define COLLECT_TRAJ_H

#include <string> 
#include "folder.h"
#include "format.h"

using namespace std;
namespace KMC
{
    class CollTraj
    {
        public:
            Folder folder;
	    Format format;
        public:
	    void colltraj(string root_dir,int iter_num, int conf_num);
    };
}
#endif  // COLLECT_TRAJ_H
