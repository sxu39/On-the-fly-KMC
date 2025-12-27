#ifdef MINIMIZE_CLASS
// clang-format off
MinimizeStyle(quickmin,minQuickMin);
// clang-format on
#else

#ifndef KMC_MIN_QUICKMIN_H
#define KMC_MIN_QUICKMIN_H

#include "min.h"

namespace KMC
{
    class minQuickMin : public min
    {
        public:
            minQuickMin();
            minQuickMin(prec_t, prec_t, prec_t dt=0.001, prec_t last_negative=0);
            ~minQuickMin() override;
            void init(cell *, potential *, prec_t, prec_t);
            void set_param(const std::vector<size_t> &, const std::vector<prec_t> &) final;
            void setup_style() final;
            void reset_vectors() final;
	    index_t iterate(size_t, size_t, const std::string &) final;
            void evolve(cell *, size_t) final;
        private:
            prec_t dt;
            bigint_t last_negative;
    };
}

#endif
#endif
