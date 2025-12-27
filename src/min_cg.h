#ifdef MINIMIZE_CLASS
// clang-format off
MinimizeStyle(cg, minCG);
// clang-format on
#else

#ifndef KMC_MIN_CG_H
#define KMC_MIN_CG_H

#include "min_linesearch.h"

namespace KMC
{
    class minCG : public minLineSearch
    {
        private:
            size_t nlimit; // max # of CG iterations before restarting
        public:
            minCG();
            minCG(prec_t tol_ene, prec_t tol_force);
            void init(cell *, potential *, size_t);
            void set_param(const std::vector<size_t> &, const std::vector<prec_t> &);
	    index_t iterate(size_t, size_t, const std::string &) override;
    };
}

#endif
#endif
