#ifndef KMC_IN_OUT_H
#define KMC_IN_OUT_H

#include "cell.h"
#include "common.h"
#include <iostream>
#include <fstream>

namespace KMC
{
    class in_out
    {
        public:
            std::string header;
            prec_t scaling;
            std::vector<std::string> elements;
        public:
            std::ofstream out_log;
            std::ofstream out_phy;
            size_t quan_wid;
        public:
            in_out();
            explicit in_out(const std::string &, const std::string &, size_t);
            ~in_out();
            void set_log(const std::string &);
            void set_phy(const std::string &);
            const cell read_stru(const std::string &);
            void print_stru(const cell &, const std::string &stru_file="POSCAR", int precision=10) const;
    };
}

#endif
