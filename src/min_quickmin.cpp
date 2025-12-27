#include "min_quickmin.h"
#include <iomanip>

using namespace KMC;
using namespace std;

// const prec_t EPS_ENERGY = 1.0e-8; // EPS_ENERGY = minimum normalization for energy tolerance
const size_t DELAYSTEP = 5;

KMC::minQuickMin::minQuickMin():min(), dt(0.001), last_negative(0){}

KMC::minQuickMin::minQuickMin(prec_t tol_ene, prec_t tol_force, prec_t dt, prec_t last_negative):
min(tol_ene, tol_force), dt(dt), last_negative(last_negative){}

KMC::minQuickMin::~minQuickMin(){}

void KMC::minQuickMin::init(cell *cell, potential *pot, prec_t dt, prec_t last_negative)
{
    min::init(pot);
    set_system(cell);
    this->dt = dt;
    this->last_negative = last_negative;
}

void KMC::minQuickMin::set_param(const vector<size_t> &iparams, const vector<prec_t> &fparams)
{
    if (fparams.size() > 0)
    {
        dt = fparams[0];
        if (fparams.size() > 1)
            last_negative = fparams[1];
    }
}

void KMC::minQuickMin::setup_style()
{
    system->velocities.set_coords();
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void KMC::minQuickMin::reset_vectors(){}

/* ----------------------------------------------------------------------
   minimization via QuickMin damped dynamics
------------------------------------------------------------------------- */

index_t KMC::minQuickMin::iterate(size_t maxiter, size_t interval, const std::string &header)
{
    size_t width = io_p->quan_wid;
    bigint_t ntimestep;
    prec_t vmax, vdotf, fdotf, scale;
    prec_t dtv, dtf, dtfm;

    coord &v = system->velocities;
    coord &f = system->forces;
    const atom_vec<bool_t> &m = system->get_moves();
    size_t atom_num = system->get_atom_num();

    for (int iter = 0; iter < maxiter; iter++)
    {
        // ntimestep = ++update->ntimestep;
        // niter++;
        ntimestep = iter;

        // zero velocity if anti-parallel to force
        // else project velocity in direction of force

        vdotf = 0.0;
        for (int i = 0; i < atom_num; i++)
            vdotf += v[i][0]*f[i][0]*m[i][0] + v[i][1]*f[i][1]*m[i][1] + v[i][2]*f[i][2]*m[i][2];

        if (vdotf < 0.0)
        {
            last_negative = ntimestep;
            for (int i = 0; i < atom_num; i++)
                v[i][0] = v[i][1] = v[i][2] = 0.0;
        }
        else
        {
            fdotf = 0.0;
            for (int i = 0; i < atom_num; i++)
                fdotf += f[i][0]*f[i][0]*m[i][0] + f[i][1]*f[i][1]*m[i][1] + f[i][2]*f[i][2]*m[i][2];

            if (fdotf == 0.0) scale = 0.0;
            else scale = vdotf/fdotf;
            for (int i = 0; i < atom_num; i++)
            {
                v[i][0] = scale * f[i][0] * m[i][0];
                v[i][1] = scale * f[i][1] * m[i][1];
                v[i][2] = scale * f[i][2] * m[i][2];
            }
        }

        // limit timestep so no particle moves further than dmax

        const vector<prec_t> &mass = system->get_mass();
        const vector<index_t> &type = system->get_atype();

        dtv = dt;

        for (int i = 0; i < atom_num; i++) {
            vmax = max(max(fabs(v[i][0]), fabs(v[i][1])), fabs(v[i][2]));
            if (dtv*vmax > dist_max) dtv = dist_max/vmax;
        }

        dtf = dtv * ftm2v;

        // Euler integration step

        coord &x = system->positions;

        for (int i = 0; i < atom_num; i++) {
            dtfm = dtf / mass[type[i]];
            x[i][0] += dtv * v[i][0];
            x[i][1] += dtv * v[i][1];
            x[i][2] += dtv * v[i][2];
            v[i][0] += dtfm * f[i][0] * m[i][0];
            v[i][1] += dtfm * f[i][1] * m[i][1];
            v[i][2] += dtfm * f[i][2] * m[i][2];
        }

        prev_ene = curr_ene;
        curr_ene = energy_force();
        if ((iter+1) % 1 == 0)
            io_p->out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
	    if (tol_ene > 0.0 && ntimestep-last_negative > DELAYSTEP)
                if (fabs(curr_ene-prev_ene) < tol_ene)
                    io_p->out_phy << "Final Energy  " << curr_ene << endl;
        neval++;

        // output structure
        if (interval)
            if ((iter+1) % interval == 0)
                io_p->print_stru(*system, header+to_string(iter+1));

        // energy tolerance criterion
        // only check after DELAYSTEP elapsed since velocties reset to 0
        // sync across replicas if running multi-replica minimization

        if (tol_ene > 0.0 && ntimestep-last_negative > DELAYSTEP)
            if (fabs(curr_ene-prev_ene) < tol_ene)
            {
                //io_p->out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
                return ETOL; // only return ETOL when the criterion is satisfied by all images
            }

        // force tolerance criterion
        // sync across replicas if running multi-replica minimization

        fdotf = 0.0;
        if (tol_force > 0.0)
        {
            if (force_style == MAX) fdotf = f.norm_max(m);       // max force norm
            else if (force_style == INF) fdotf = f.norm_inf(m);  // inf force norm
            else if (force_style == TWO) fdotf = f.norm_sqr(m);  // Euclidean force 2-norm
            else throw "Illegal min_modify command.";
            if (fdotf < tol_force*tol_force)
            {
                io_p->out_phy << setprecision(9) << setw(width) << iter + 1 << setw(width) << curr_ene << "\n";
                return FTOL; // only return FTOL when the criterion is satisfied by all images
            }
        }
    }
    return MAXITER;
}

void KMC::minQuickMin::evolve(cell *system, size_t ntimestep)
{
    // zero velocity if anti-parallel to force
    // else project velocity in direction of force

    coord &v = system->velocities;
    coord &f = system->forces;
    const atom_vec<bool_t> &m = system->get_moves();
    size_t atom_num = system->get_atom_num();

    prec_t vmax, vdotf, fdotf, scale;
    prec_t dtv, dtf, dtfm;

    vdotf = 0.0;
    for (int i = 0; i < atom_num; i++)
        vdotf += v[i][0]*f[i][0]*m[i][0] + v[i][1]*f[i][1]*m[i][1] + v[i][2]*f[i][2]*m[i][2];

    if (vdotf < 0.0)
    {
        last_negative = ntimestep;
        for (int i = 0; i < atom_num; i++)
            v[i][0] = v[i][1] = v[i][2] = 0.0;
    }
    else
    {
        fdotf = 0.0;
        for (int i = 0; i < atom_num; i++)
            fdotf += f[i][0]*f[i][0]*m[i][0] + f[i][1]*f[i][1]*m[i][1] + f[i][2]*f[i][2]*m[i][2];

        if (fdotf == 0.0) scale = 0.0;
        else scale = vdotf/fdotf;
        for (int i = 0; i < atom_num; i++)
        {
            v[i][0] = scale * f[i][0] * m[i][0];
            v[i][1] = scale * f[i][1] * m[i][1];
            v[i][2] = scale * f[i][2] * m[i][2];
        }
    }

    // limit timestep so no particle moves further than dmax

    const vector<prec_t> &mass = system->get_mass();
    const vector<index_t> &type = system->get_atype();

    dtv = dt;

    for (int i = 0; i < atom_num; i++) {
        vmax = max(max(fabs(v[i][0]), fabs(v[i][1])), fabs(v[i][2]));
        if (dtv*vmax > dist_max) dtv = dist_max/vmax;
    }

    dtf = dtv * ftm2v;

    // Euler integration step

    coord &x = system->positions;

    for (int i = 0; i < atom_num; i++)
    {
        dtfm = dtf / mass[type[i]];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        v[i][0] += dtfm * f[i][0] * m[i][0];
        v[i][1] += dtfm * f[i][1] * m[i][1];
        v[i][2] += dtfm * f[i][2] * m[i][2];
    }
}
