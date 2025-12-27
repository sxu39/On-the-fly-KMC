#ifndef CONF_SCREEN_H
#define CONF_SCREEN_H

#include "cell.h"
#include "in_out.h"
#include <string>

using namespace std;

namespace KMC
{
    class ConfScreen
    {
        public:
            in_out io_obj;
            cell opt_cell;
        public:
            void conf_screen(string iniconfrlx_dir, string confcand_folder_dir, string conf_folder_dir, double dist_threshold);
    };
}
#endif
