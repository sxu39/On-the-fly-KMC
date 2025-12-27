#include "min.h"
using namespace KMC;
using namespace std;

KMC::min::min():searchflag(0), io_p(nullptr), system(nullptr), pot(nullptr), dist_max(0.2), tol_ene(0.0), tol_force(0.0), prev_ene(0.0),
curr_ene(0.0), force_style(TWO), alpha(0.0), max_eval(INFINITY){}

KMC::min::min(prec_t tol_ene, prec_t tol_force):searchflag(0), io_p(nullptr), system(nullptr), pot(nullptr), dist_max(0.2), tol_ene(tol_ene), 
tol_force(tol_force), prev_ene(0.0), curr_ene(0.0), force_style(TWO), alpha(0.0), max_eval(INFINITY){}

void KMC::min::init(potential *pot)
{
    this->pot = pot;
}

void KMC::min::set_in_out(in_out *io_p)
{
    this->io_p = io_p;
}

void KMC::min::set_system(cell *cell)
{
    system = cell;
    niter = 0;
    neval = 0;
    curr_ene = energy_force();
    setup_style();
    reset_vectors();
}

KMC::min::~min(){}

prec_t KMC::min::energy_force()
{
    prec_t energy = pot->infer(system);
    return energy;
}

char *KMC::min::stopstrings(index_t n)
{
  const char *strings[] = {"max iterations",
                           "max force evaluations",
                           "energy tolerance",
                           "force tolerance",
                           "search direction is not downhill",
                           "linearsearch alpha is zero",
                           "forces are zero",
                           "quadratic factors are zero"};
  return (char *) strings[n];
}
