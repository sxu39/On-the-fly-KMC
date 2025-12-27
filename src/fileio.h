#ifndef FILEIO_H
#define FILEIO_H

#include "para.h"
#include <string>
#include <iostream>

using namespace std;    

class FileIO 
{
    public:
        Para para;
    public:
        Para read_status(string file_name);
	Para read_model_devi(string file_name, double aveforce_limit);
};
#endif
