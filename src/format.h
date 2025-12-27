#ifndef FORMAT_H
#define FORMAT_H

#include <string> 
#include "in_out.h"
#include "cell.h"

using namespace std;

namespace KMC
{
    class Format
    {
	public:
            in_out io_obj;
            cell opt_cell;
	public:
	    void pos2dump(string poscar_path, string dump_path);
    };
}

#endif  // FORMAT_H
