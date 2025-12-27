#ifndef KMC_MIN_LSRCH_H
#define KMC_MIN_LSRCH_H

#include "min.h"

namespace KMC
{
    class minLineSearch : public min
    {
        public:
            minLineSearch();
            minLineSearch(prec_t, prec_t);
            ~minLineSearch() override;
            void init(potential *) override;
            void setup_style() override;
            void reset_vectors() override;

        protected:
            // vectors needed by linesearch minimizers
            // allocated and stored by fix_minimize
            // x,f are stored by parent or Atom class or Pair class

            std::vector<prec_t> x0;    // coords at start of linesearch
            prec_t *g;     // old gradient vector
            prec_t *h;     // search direction vector

            size_t nvec;

            prec_t *xvec;   // variables for atomic dof, as 1d vector
            prec_t *fvec;   // force vector for atomic dof, as 1d vector
            bool_t *mvec;   // variables for atom mobility, as 1d vector

            typedef index_t (minLineSearch::*FnPtr)(prec_t, prec_t &);
            FnPtr linemin;
            index_t linemin_backtrack(prec_t, prec_t &);
            index_t linemin_quadratic(prec_t, prec_t &);
            index_t linemin_forcezero(prec_t, prec_t &);
            void evolve(cell *, size_t) final {}

            prec_t alpha_step(prec_t);
            prec_t compute_dir_deriv(prec_t &);
    };
}

#endif
