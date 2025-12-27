#ifndef KMC_OPT_H
#define KMC_OPT_H

#include "min.h"
#include "in_out.h"
#include "timer.h"

namespace KMC
{
    class opt
    {
        private:
            cell *system;
            min *minimize;
            MinimizeCreatorMap *minimize_map;
            potential *pot;
            timer *timer_p;
            prec_t tol_ene;
            prec_t tol_force;
        public:
            void set_minimize(const std::string &);
            void set_potential(const std::string &);
        public:
            opt();
            opt(prec_t, prec_t);
            ~opt();
            void init(const std::string &, const std::string &, in_out *, timer *, const std::vector<size_t> &iparams=std::vector<size_t>(), 
            const std::vector<prec_t> &fparams=std::vector<prec_t>());
            void set_system(cell *);
            void run(size_t, size_t interval=0, const std::string &label="");
    };
}

#endif
