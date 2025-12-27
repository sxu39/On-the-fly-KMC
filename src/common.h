#ifndef KMC_COMMON_H
#define KMC_COMMON_H

#include "type.h"
#include <unordered_map>
#include <string>

namespace KMC
{
    const prec_t kB = 1.38e-23;
    const prec_t ftm2v = 1.0 / 1.0364269e-4;
    const std::unordered_map<std::string, prec_t> ele_mass{{"H", 1.008}, {"O", 16.00}, {"Li", 6.941}, {"C", 12.01}, {"Cu", 63.55}, {"Pt", 195.08}, {"P", 30.974},{"F",18.998}};
    const double pai = 3.141592653;
    const double convert_factor = 9648.68; 
    //const double T = 300.0;
    const double NA = 6.02e23;
    const double dt = 0.001;//time step unit:ps
}

#endif
