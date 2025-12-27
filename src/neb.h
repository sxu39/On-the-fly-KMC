#ifndef KMC_NEB_H
#define KMC_NEB_H

#include "common.h"
#include "min.h"
#include "timer.h"
#include "utils.h"

namespace KMC {
    struct image
    {
        prec_t nlen = 0.0;
        prec_t plen = 0.0;
        prec_t tlen = 0.0;
        prec_t gradlen = 0.0;
        prec_t gradnextlen = 0.0;
        prec_t dotpath = 0.0;
        prec_t dottangrad = 0.0;
        prec_t dotgrad = 0.0;
        coord tangent;
        coord springF;
        prec_t dot = 0.0;
        prec_t dotSpringTangent = 0.0;
    };

    class neb
    {
        private:
            cell *initial_system;
            cell *final_system;
            std::vector<cell> inner_system;
            std::vector<prec_t> curr_ene;
            std::vector<prec_t> prev_ene;
            std::vector<image> sys_info;
            index_t n_climber;
            timer *timer_p;
            min *minimize;
            prec_t dt;
            bigint_t last_negative;
            potential *pot;
            size_t n_image;
            prec_t k_spring;
            prec_t k_spring_perp;
            prec_t time_step;
            prec_t tol_ene;
            prec_t tol_force;
            size_t n1_steps;
            size_t n2_steps;
        private:
            //void interpolate();
            void post_force_calculate(cell *, cell *, cell *, index_t, cell_type);
            void post_force_assign(cell *, index_t, cell_type);
            void post_force();
	    index_t iterate(size_t, size_t, const std::string &);
        public:
            neb();
            ~neb();
            inline void set_timer(timer *);
            inline void set_in_out(in_out *);
	    void interpolate();
            void set_system(cell *, cell *);
            void set_potential(const std::string &);
	    void run(size_t interval=0, const std::string &label="");
	    //void out_stru(const std::string &header="POSCAR_") const;
	    void out_stru(std::string temp_pos_path) const;
    };
}

void KMC::neb::set_timer(timer *Timer)
{
    timer_p = Timer;
}

void KMC::neb::set_in_out(in_out *IO)
{
    minimize->set_in_out(IO);
}

#endif
