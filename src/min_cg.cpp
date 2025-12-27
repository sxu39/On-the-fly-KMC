#include "min_cg.h"
#include <cmath>
#include <iomanip>

using namespace KMC;
using namespace std;

// EPS_ENERGY = minimum normalization for energy tolerance

// const prec_t EPS_ENERGY = 1.0e-8;

/* ---------------------------------------------------------------------- */

KMC::minCG::minCG():minLineSearch(), nlimit(10){}

KMC::minCG::minCG(prec_t tol_ene, prec_t tol_force):minLineSearch(tol_ene, tol_force), nlimit(10){}

void KMC::minCG::init(cell *cell, potential *pot, size_t nlimit)
{
    minLineSearch::init(pot);
    set_system(cell);
    this->nlimit = nlimit;
}

void KMC::minCG::set_param(const vector<size_t> &iparams, const vector<prec_t> &fparams)
{
    if (iparams.size() > 0)
        nlimit = iparams[0];
}

/* ----------------------------------------------------------------------
   minimization via conjugate gradient iterations
------------------------------------------------------------------------- */

index_t KMC::minCG::iterate(size_t maxiter, size_t interval, const string &header)
{
    size_t width = io_p->quan_wid;
    index_t i, fail;
    prec_t beta, gg, dot[2], fdotf;
    const bool_t *m = mvec;

    // initialize working vectors

    for (i = 0; i < nvec; i++) h[i] = g[i] = fvec[i] * m[i];

    gg = system->forces.norm_sqr();

    for (int iter = 0; iter < maxiter; iter++)
    {

        // ntimestep = ++update->ntimestep;
        // niter++;

        // line minimization along direction h from current atom->x

        prev_ene = curr_ene;
        fail = (this->*linemin)(curr_ene, alpha);
        if (fail) return fail;

        // function evaluation criterion

        if (neval >= max_eval) return MAXEVAL;
        
        if ((iter+1) % 1 == 0)
            io_p->out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
            if (fabs(curr_ene-prev_ene) < tol_ene)
                io_p->out_phy << "Final Energy  " << curr_ene << endl;
	else if (iter == maxiter - 1)
            io_p->out_phy << "Final Energy  " << curr_ene << endl;

        // output structure
        if (interval)
            if (iter % interval == 0)
                io_p->print_stru(*system, header+to_string(iter+1));

        // energy tolerance criterion

        if (fabs(curr_ene-prev_ene) < tol_ene)
        {
            //out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
            return ETOL;
        }

        // force tolerance criterion

        dot[0] = dot[1] = 0.0;
        for (i = 0; i < nvec; i++) {
            dot[0] += fvec[i] * fvec[i] * m[i];
            dot[1] += fvec[i] * g[i];
        }

        fdotf = 0.0;
        if (tol_force > 0.0)
        {
            if (force_style == MAX) fdotf = system->forces.norm_max(system->get_moves());       // max force norm
            else if (force_style == INF) fdotf = system->forces.norm_inf(system->get_moves());  // inf force norm
            else if (force_style == TWO) fdotf = system->forces.norm_sqr(system->get_moves());  // Euclidean force 2-norm
            else throw "Illegal min_modify command.";
            if (fdotf < tol_force*tol_force)
            {    
                io_p->out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
                return FTOL; // only return FTOL when the criterion is satisfied by all images
            }
        }

        // update new search direction h from new f = -Grad(x) and old g
        // this is Polak-Ribieri formulation
        // beta = dotall[0]/gg would be Fletcher-Reeves
        // reinitialize CG every ndof iterations by setting beta = 0.0

        beta = max(0.0,(dot[0] - dot[1])/gg);
        if ((niter+1) % nlimit == 0) beta = 0.0;
        gg = dot[0];

        for (i = 0; i < nvec; i++) {
            g[i] = fvec[i] * m[i];
            h[i] = g[i] + beta * h[i];
        }
        // reinitialize CG if new search direction h is not downhill

        dot[0] = 0.0;
        for (i = 0; i < nvec; i++) dot[0] += g[i] * h[i];

        if (dot[0] <= 0.0) {
            for (i = 0; i < nvec; i++) h[i] = g[i];
        }
    }
    return MAXITER;
}
