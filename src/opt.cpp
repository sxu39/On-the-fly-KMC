#include "opt.h"
#include "opt_minimize.h"
#include <iomanip>
using namespace KMC;
using namespace std;

KMC::opt::opt():system(nullptr), minimize(nullptr), pot(nullptr), timer_p(nullptr), tol_ene(1e-6), tol_force(0)
{
    minimize_map = new MinimizeCreatorMap;
    #define MINIMIZE_CLASS
    #define MinimizeStyle(key, Class) (*minimize_map)[#key] = &minimize_creator<Class>;
    #include "opt_minimize.h"
    #undef MinimizeStyle
    #undef MINIMIZE_CLASS
}

KMC::opt::opt(prec_t tol_ene, prec_t tol_force):system(nullptr), minimize(nullptr), pot(nullptr), timer_p(nullptr), tol_ene(tol_ene), tol_force(tol_force)
{
    minimize_map = new MinimizeCreatorMap;
    #define MINIMIZE_CLASS
    #define MinimizeStyle(key, Class) (*minimize_map)[#key] = &minimize_creator<Class>;
    #include "opt_minimize.h"
    #undef MinimizeStyle
    #undef MINIMIZE_CLASS
}

KMC::opt::~opt()
{
    if (minimize) delete minimize;
    if (pot) delete pot;
    delete minimize_map;
}

void KMC::opt::init(const string &pot_name, const string &min_style, in_out *io_obj, timer *timer_obj, const vector<size_t> &iparams, const vector<prec_t> &fparams)
{
    set_potential(pot_name);
    set_minimize(min_style);
    minimize->init(pot);
    minimize->set_in_out(io_obj);
    timer_p = timer_obj;
    minimize->set_param(iparams, fparams);
}

void KMC::opt::set_system(cell *system)
{
    this->system = system;
    minimize->set_system(system);
}

void KMC::opt::set_minimize(const string &min_style)
{
    try
    {
        if (minimize) delete minimize;
        if (minimize_map->find(min_style) != minimize_map->end())
        {
            MinimizeCreator &minimize_creator = (*minimize_map)[min_style];
            minimize = minimize_creator(tol_ene, tol_force);
        }
        else throw "Illegal minimize style.";
    }
    catch (const char *msg)
    {
        cerr << msg << "\n";
    }
}

void KMC::opt::set_potential(const string &pot_name)
{
    pot = new ml_potential(pot_name);
}

void KMC::opt::run(size_t maxiter, size_t interval, const string &label)
{
    in_out *io_p = minimize->io_p;
    string header;
    //string header = "POSCAR_OPT_";
    //if (label!="") header += label + "_";
    if (label!="") header = label + "_";
    // minimizer iterations

    io_p->out_log << timer_p->current_time() << "single structure minimization" << std::endl;
    io_p->out_phy << "single structure minimization" << std::endl;
    io_p->out_phy << setw(io_p->quan_wid) << "step" << setw(io_p->quan_wid) << "Energy" << std::endl;
    size_t stop_condition = minimize->iterate(maxiter, interval, header);
    char *stopstr = min::stopstrings(stop_condition);

    if (stop_condition != MAXITER){
        io_p->out_log << timer_p->current_time() << stopstr << "\n";
    }
    io_p->out_log << timer_p->current_time() << "OPT END" << std::endl;
}
