#include <iostream>
#include <fstream>
#include <sstream>
#include "fileio.h"

Para FileIO::read_status(string file_name) 
{
    ifstream ifs(file_name);
    if (!ifs.is_open()) 
    {
        int num_yes = 0;
        para.yes_count = num_yes;
        ofstream ofs(file_name);
        if (!ofs.is_open()) 
        {
            cerr << "Can't create status file" << std::endl;
        }
        ofs.close();
    } 
    else 
    {
        string line;
        int num_yes = 0;
        while (getline(ifs, line)) {
            if (!line.empty()) {
                size_t pos = 0;
                while ((pos = line.find("YES", pos)) != string::npos) {
                    num_yes += 1;
                    pos += 3;
                }
            }
        }
        para.yes_count = num_yes;
    }
    ifs.close();
    return para;
}

Para FileIO::read_model_devi(string file_name, double aveforce_limit)
{
    ifstream ifs(file_name);
    if (!ifs.is_open())
    {
        double devi = aveforce_limit;
        para.upper_devi = devi;
        ofstream ofs(file_name);
        if (!ofs.is_open())
        {
            cerr << "Can't create status file" << endl;
        }
        ofs.close();
    }
    else
    {
        string line;
        string lastLine;
        double devi;
        while (getline(ifs, line))
        {
            lastLine = line;
        }
        istringstream iss(lastLine);
        int value1;
        double value2;
        if (iss >> value1)
        {
            if (iss >> value2)
            {
                devi = value2;
            }
        }
        para.upper_devi = devi;
    }
    ifs.close();
    return para;
}
