#ifndef KMC_MIN_H
#define KMC_MIN_H

#include "cell.h"
#include "common.h"
#include "potential.h"
#include "in_out.h"
#include <cmath>
#include <map>

namespace KMC
{
    enum {MAXITER, MAXEVAL, ETOL, FTOL, DOWNHILL, ZEROALPHA, ZEROFORCE, ZEROQUAD};
    class min
    {
        friend class neb;
        friend class opt;
        protected:
            index_t searchflag; // 0 if damped dynamics, 1 if sub-cycles on local search
            in_out *io_p;
            cell *system;
            potential *pot;
            prec_t dist_max; // max distance to move any atom in one step
            prec_t tol_ene;
            prec_t tol_force;
            size_t niter;
            size_t neval;
            prec_t prev_ene;
            prec_t curr_ene;
            norm_style force_style;
            prec_t alpha;
        public:
            size_t max_eval;
        private:
            static char *stopstrings(index_t n);
        protected:
            prec_t energy_force();
            virtual index_t iterate(size_t, size_t, const std::string &) = 0;
        public:
            min();
            min(prec_t, prec_t);
            virtual ~min();
            virtual void init(potential *);
            void set_in_out(in_out *);
            void set_system(cell *);
            virtual void evolve(cell *, size_t) = 0;
            virtual void set_param(const std::vector<size_t> &, const std::vector<prec_t> &) = 0;
            virtual void setup_style() = 0;
            virtual void reset_vectors() = 0;
    };

    typedef min *(*MinimizeCreator)(prec_t, prec_t);
    typedef std::map<std::string, MinimizeCreator> MinimizeCreatorMap;
}

namespace KMC {
    template <class T> static min *minimize_creator(prec_t tol_ene, prec_t tol_force)
    {
        return new T(tol_ene, tol_force);
    }
}

#endif
