#include "neb.h"
#include "min_quickmin.h"
#include <iomanip>
using namespace KMC;
using namespace std;

const size_t DELAYSTEP = 5;
const prec_t MY_PI = acos(-1);

KMC::neb::neb():initial_system(nullptr), final_system(nullptr), pot(nullptr), timer_p(nullptr), n_image(5), k_spring(10.0), k_spring_perp(0.0), 
tol_ene(0.005), tol_force(0.001), n1_steps(100), n2_steps(100), dt(0.001), last_negative(0.0)
{
    minimize = new minQuickMin(tol_ene, tol_force, dt, last_negative);
    inner_system.resize(n_image-2);
    prev_ene.resize(n_image);
    curr_ene.resize(n_image);
    sys_info.resize(n_image);
}

void KMC::neb::set_system(cell *initial, cell *final)
{
    initial_system = initial;
    final_system = final;
    sys_info[0].tangent.init_coords(initial_system->get_atom_num()*3);
    sys_info[n_image-1].tangent.init_coords(final_system->get_atom_num()*3);
    sys_info[0].springF.init_coords(initial_system->get_atom_num()*3);
    sys_info[n_image-1].springF.init_coords(final_system->get_atom_num()*3);
}

void KMC::neb::set_potential(const string &pot_name)
{
    pot = new ml_potential(pot_name);
}

KMC::neb::~neb()
{
    delete minimize;
    if (pot) delete pot;
}

void KMC::neb::interpolate()
{
    prec_t delx, dely, delz;
    size_t n_atom = initial_system->get_atom_num();

    vector<prec_t> del_position(3*n_atom);
    for (int i = 0 ; i < n_atom ; i++)
    {
        delx = final_system->positions[i][0] - initial_system->positions[i][0];
        dely = final_system->positions[i][1] - initial_system->positions[i][1];
        delz = final_system->positions[i][2] - initial_system->positions[i][2];
        initial_system->minimum_image(delx, dely, delz);
        del_position[3*i] = delx;
        del_position[3*i+1] = dely;
        del_position[3*i+2] = delz;
    }

    for (int i_image = 1 ; i_image < n_image - 1 ; i_image++)
    {
        prec_t fraction = i_image / (n_image - 1.0);
        inner_system[i_image-1].init(*initial_system);
        sys_info[i_image].tangent.init_coords(initial_system->get_atom_num()*3);
        sys_info[i_image].springF.init_coords(initial_system->get_atom_num()*3);
        for (int i = 0 ; i < n_atom ; i++)
            for (int j = 0 ; j < 3 ; j++)
                inner_system[i_image-1].positions[i][j] += fraction * del_position[3*i+j];
    }
}

void KMC::neb::run(size_t interval, const string &label)
{
    in_out *io_p = minimize->io_p;
    string header = "POSCAR_NEB_";
    if (label!="") header += label + "_";
    //interpolate();
    size_t stop_condition;
    char *stopstr;

    io_p->out_log << timer_p->current_time() << "NEB minimization" << std::endl;
    io_p->out_phy << "NEB minimization" << std::endl;
    io_p->out_phy << setw(io_p->quan_wid) << "step";
    for (int i = 0 ; i < n_image ; i++)
        io_p->out_phy << setw(io_p->quan_wid) << "E_image_"+to_string(i);
    io_p->out_phy << std::endl;
    stop_condition = iterate(n1_steps, interval, header);
    stopstr = min::stopstrings(stop_condition);
    if (stop_condition != MAXITER)
    {
        io_p->out_log << timer_p->current_time() << stopstr << "\n";
    }

    header = "POSCAR_CINEB_";
    if (label!="") header += label + "_";

    prec_t max_ene = curr_ene[0];
    for (int i = 1 ; i < n_image ; i++)
        if (max_ene < curr_ene[i])
        {
            max_ene = curr_ene[i];
            n_climber = i;
        }

    io_p->out_log << timer_p->current_time() << "CI-NEB minimization" << endl;
    io_p->out_phy << "CI-NEB minimization" << endl;
    io_p->out_phy << setw(io_p->quan_wid) << "step";
    for (int i = 0 ; i < n_image ; i++)
        io_p->out_phy << setw(io_p->quan_wid) << "E_image_"+to_string(i);
    io_p->out_phy << endl;
    stop_condition = iterate(n2_steps, interval, header);
    stopstr = min::stopstrings(stop_condition);
    if (stop_condition != MAXITER)
    {
        io_p->out_log << timer_p->current_time() << stopstr << "\n";
    }
    io_p->out_log << timer_p->current_time() << "NEB END" << std::endl;
}

index_t KMC::neb::iterate(size_t maxiter, size_t interval, const string &header)
{
    ofstream &out_phy = minimize->io_p->out_phy;
    size_t width = minimize->io_p->quan_wid;
    vector<prec_t> inner_enes;
    curr_ene[0] = pot->infer(initial_system);
    inner_enes = pot->infer(inner_system);
    for (int i = 1 ; i < n_image - 1 ; i++)
        curr_ene[i] = inner_enes[i-1];
    curr_ene[n_image-1] = pot->infer(final_system);
    post_force();

    minimize->neval = 0;

    for (int iter = 0; iter < maxiter; iter++)
    {
        minimize->evolve(initial_system, iter);
        for (auto &system : inner_system)
            minimize->evolve(&system, iter);
        minimize->evolve(final_system, iter);

        for (int i = 0 ; i < n_image ; i++)
            prev_ene[i] = curr_ene[i];

        curr_ene[0] = pot->infer(initial_system);
        inner_enes = pot->infer(inner_system);
        for (int i = 1 ; i < n_image - 1 ; i++)
            curr_ene[i] = inner_enes[i-1];
        curr_ene[n_image-1] = pot->infer(final_system);
        post_force();

        if ((iter+1) % 1 == 0)
        {
            out_phy << setw(width) << iter + 1 << setprecision(6);
            for (auto &curr : curr_ene) 
                out_phy << setw(width) << curr;
            out_phy << "\n";
        }
        minimize->neval++;

	// output structure
        //if (interval)
        //    if ((iter+1) % interval == 0)
        //        out_stru(header+to_string(iter+1)+"_");

        if (minimize->neval > maxiter)
            return MAXEVAL;

        // energy tolerance criterion
        // only check after DELAYSTEP elapsed since velocties reset to 0
        // sync across replicas if running multi-replica minimization

        if (tol_ene > 0.0 && iter-last_negative > DELAYSTEP)
     	{
            if (max_vec(fabs(curr_ene-prev_ene)) < tol_ene)
            {
                //out_phy << setw(width) << iter + 1 << setprecision(6);
                ofstream out_ene("neb_ene.txt");	
	        out_phy << setw(width) << "Final Energy " << setprecision(6);
	        int num = 0;
                for (auto &curr : curr_ene)
        	{	
		    out_ene << setprecision(width) << num << "   " << curr << endl;
                    out_phy << setw(width) << curr;
		    num++;
	        }
                out_phy << "\n";
        	/*
        	for (auto &curr : curr_ene)
        	{
                    out_ene << curr << endl;
	        }
	        */
                return ETOL; // only return ETOL when the criterion is satisfied by all images
            }
    	}
	if ((iter + 1) == maxiter)
        {
            ofstream out_ene("neb_ene.txt");
            out_phy << setw(width) << "Final Energy " << setprecision(6);
            int num = 0;
            for (auto &curr : curr_ene)
            {
                out_ene << setprecision(width) << num << "   " << curr << endl;
                out_phy << setw(width) << curr;
                num++;
            }
            out_phy << "\n";
        }

        // force tolerance criterion
        // sync across replicas if running multi-replica minimization

        vector<prec_t> fdotf(n_image);
        if (tol_force > 0.0)
        {
            if (minimize->force_style == MAX)
            {
                fdotf[0] = initial_system->forces.norm_max(initial_system->moves);
                for (int i = 1 ; i < n_image - 1 ; ++i)
                    fdotf[i] = inner_system[i-1].forces.norm_max(inner_system[i-1].moves);
                fdotf[n_image-1] = final_system->forces.norm_max(final_system->moves);
            }  // max force norm
            else if (minimize->force_style == INF)
            {
                fdotf[0] = initial_system->forces.norm_inf(initial_system->moves);
                for (int i = 1 ; i < n_image - 1 ; ++i)
                    fdotf[i] = inner_system[i-1].forces.norm_inf(inner_system[i-1].moves);
                fdotf[n_image-1] = final_system->forces.norm_inf(final_system->moves);
            }  // inf force norm
            else if (minimize->force_style == TWO)
            {
                fdotf[0] = initial_system->forces.norm_sqr(initial_system->moves);
                for (int i = 1 ; i < n_image - 1 ; ++i)
                    fdotf[i] = inner_system[i-1].forces.norm_sqr(inner_system[i-1].moves);
                fdotf[n_image-1] = final_system->forces.norm_sqr(final_system->moves);
            }  // Euclidean force 2-norm
            else throw "Illegal min_modify command.";
            if (max_vec(fdotf) < tol_force*tol_force)
            {
                out_phy << setw(width) << iter + 1 << setprecision(9);
                for (auto &curr : curr_ene) 
                    out_phy << setw(width) << curr;
                out_phy << "\n";
                return FTOL; // only return FTOL when the criterion is satisfied by all images
            }
        }
    }
    return MAXITER;
}

void KMC::neb::post_force()
{
    post_force_calculate(initial_system, nullptr, &inner_system[0], 0, INITIAL);
    post_force_calculate(&inner_system[0], initial_system, &inner_system[1], 1, INNER);
    for (int i = 1 ; i < n_image - 3 ; ++i)
        post_force_calculate(&inner_system[i], &inner_system[i-1], &inner_system[i+1], i+1, INNER);
    post_force_calculate(&inner_system[n_image-3], &inner_system[n_image-4], final_system, n_image-2, INNER);
    post_force_calculate(final_system, &inner_system[n_image-3], nullptr, n_image-1, FINAL);

    post_force_assign(initial_system, 0, INITIAL);
    for (int i = 1 ; i < n_image - 1 ; i++)
    {
        if (i == n_climber)
            post_force_assign(&inner_system[i-1], i, CLIMBER);
        else
            post_force_assign(&inner_system[i-1], i, INNER);
    }
    post_force_assign(final_system, n_image-1, FINAL);
}

void KMC::neb::post_force_calculate(cell *system, cell *prev_system, cell *next_system, index_t sys_id, cell_type system_type)
{
    prec_t vprev, vnext, veng;
    veng = curr_ene[sys_id];

    coord &x = system->positions;
    coord &f = system->forces;
    size_t atom_num = system->get_atom_num();
    prec_t delxp, delyp, delzp, delxn, delyn, delzn;

    prec_t plen = 0.0;
    prec_t nlen = 0.0;
    prec_t tlen = 0.0;
    prec_t gradnextlen = 0.0;

    prec_t dotgrad, gradlen, dotpath, dottangrad;
    dotgrad = gradlen = dotpath = dottangrad = 0.0;

    coord &tangent = sys_info[sys_id].tangent;
    coord &springF = sys_info[sys_id].springF;

    if (system_type == FINAL)
    {
        coord &xprev = prev_system->positions;
        for (int i = 0; i < atom_num; i++)
        {
            delxp = x[i][0] - xprev[i][0];
            delyp = x[i][1] - xprev[i][1];
            delzp = x[i][2] - xprev[i][2];
            system->minimum_image(delxp, delyp, delzp);
            plen += delxp * delxp + delyp * delyp + delzp * delzp;
            dottangrad += delxp * f[i][0] + delyp * f[i][1] + delzp * f[i][2];
            gradlen += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
        }

    }
    else if (system_type == INITIAL)
    {
        coord &xnext = next_system->positions;
        coord &fnext = next_system->forces;
        for (int i = 0; i < atom_num; i++)
        {
            delxn = xnext[i][0] - x[i][0];
            delyn = xnext[i][1] - x[i][1];
            delzn = xnext[i][2] - x[i][2];
            system->minimum_image(delxn, delyn, delzn);
            nlen += delxn * delxn + delyn * delyn + delzn * delzn;
            gradnextlen += fnext[i][0] * fnext[i][0] + fnext[i][1] * fnext[i][1] + fnext[i][2] * fnext[i][2];
            dotgrad += f[i][0] * fnext[i][0] + f[i][1] * fnext[i][1] + f[i][2] * fnext[i][2];
            dottangrad += delxn * f[i][0] + delyn * f[i][1] + delzn * f[i][2];
            gradlen += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
        }
    }
    else
    {
        coord &xprev = prev_system->positions;
        coord &xnext = next_system->positions;
        coord &fnext = next_system->forces;
        vprev = curr_ene[sys_id-1];
        vnext = curr_ene[sys_id+1];
        // not the first or last replica

        double vmax = std::max(std::fabs(vnext - veng), std::fabs(vprev - veng));
        double vmin = std::min(std::fabs(vnext - veng), std::fabs(vprev - veng));

        for (int i = 0; i < atom_num; i++)
        {
            delxp = x[i][0] - xprev[i][0];
            delyp = x[i][1] - xprev[i][1];
            delzp = x[i][2] - xprev[i][2];
            system->minimum_image(delxp, delyp, delzp);
            plen += delxp * delxp + delyp * delyp + delzp * delzp;

            delxn = xnext[i][0] - x[i][0];
            delyn = xnext[i][1] - x[i][1];
            delzn = xnext[i][2] - x[i][2];
            system->minimum_image(delxn, delyn, delzn);

            if (vnext > veng && veng > vprev)
            {
                tangent[i][0] = delxn;
                tangent[i][1] = delyn;
                tangent[i][2] = delzn;
            }
            else if (vnext < veng && veng < vprev)
            {
                tangent[i][0] = delxp;
                tangent[i][1] = delyp;
                tangent[i][2] = delzp;
            }
            else
            {
                if (vnext > vprev)
                {
                    tangent[i][0] = vmax * delxn + vmin * delxp;
                    tangent[i][1] = vmax * delyn + vmin * delyp;
                    tangent[i][2] = vmax * delzn + vmin * delzp;
                }
                else if (vnext < vprev)
                {
                    tangent[i][0] = vmin * delxn + vmax * delxp;
                    tangent[i][1] = vmin * delyn + vmax * delyp;
                    tangent[i][2] = vmin * delzn + vmax * delzp;
                }
                else
                {    // vnext == vprev, e.g. for potentials that do not compute an energy
                    tangent[i][0] = delxn + delxp;
                    tangent[i][1] = delyn + delyp;
                    tangent[i][2] = delzn + delzp;
                }
            }

            nlen += delxn * delxn + delyn * delyn + delzn * delzn;
            tlen += tangent[i][0] * tangent[i][0] + tangent[i][1] * tangent[i][1] + tangent[i][2] * tangent[i][2];
            gradlen += f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
            dotpath += delxp * delxn + delyp * delyn + delzp * delzn;
            dottangrad += tangent[i][0] * f[i][0] + tangent[i][1] * f[i][1] + tangent[i][2] * f[i][2];
            gradnextlen += fnext[i][0] * fnext[i][0] + fnext[i][1] * fnext[i][1] + fnext[i][2] * fnext[i][2];
            dotgrad += f[i][0] * fnext[i][0] + f[i][1] * fnext[i][1] + f[i][2] * fnext[i][2];

            springF[i][0] = k_spring_perp * (delxn - delxp);
            springF[i][1] = k_spring_perp * (delyn - delyp);
            springF[i][2] = k_spring_perp * (delzn - delzp);
        }
    }

    sys_info[sys_id].nlen = sqrt(nlen);
    sys_info[sys_id].plen = sqrt(plen);
    sys_info[sys_id].tlen = sqrt(tlen);
    sys_info[sys_id].gradlen = sqrt(gradlen);
    sys_info[sys_id].gradnextlen = sqrt(gradnextlen);
    sys_info[sys_id].dotpath = dotpath;
    sys_info[sys_id].dottangrad = dottangrad;
    sys_info[sys_id].dotgrad = dotgrad;

    // normalization of tangent should be performed before calculation the dot product between force and tangent
    if (sys_info[sys_id].tlen > 0.0) {
        prec_t tleninv = 1.0 / sys_info[sys_id].tlen;
        for (int i = 0; i < atom_num; i++)
        {
            tangent[i][0] *= tleninv;
            tangent[i][1] *= tleninv;
            tangent[i][2] *= tleninv;
        }
    }
    
    if (system_type == INITIAL || system_type == FINAL) return;

    prec_t dot = 0.0;
    prec_t dotSpringTangent = 0.0;

    for (int i = 0; i < atom_num; i++)
    {
        dot += f[i][0] * tangent[i][0] + f[i][1] * tangent[i][1] + f[i][2] * tangent[i][2];
        dotSpringTangent += springF[i][0] * tangent[i][0] + springF[i][1] * tangent[i][1] +
            springF[i][2] * tangent[i][2];
    }
    sys_info[sys_id].dot = dot;
    sys_info[sys_id].dotSpringTangent = dotSpringTangent;
}

void KMC::neb::post_force_assign(cell *system, index_t sys_id, cell_type system_type)
{
    // normalize tangent vector

    size_t atom_num = system->get_atom_num();
    coord &tangent = sys_info[sys_id].tangent;
    coord &springF = sys_info[sys_id].springF;
    coord &f = system->forces;

    // prec_t dottangrad, dotgrad;

    // // first or last replica has no change to forces, just return

    if (system_type == INITIAL || system_type == FINAL) return;

    prec_t AngularContr;
    prec_t dotpath = sys_info[sys_id].dotpath / (sys_info[sys_id].plen * sys_info[sys_id].nlen);
    AngularContr = 0.5 * (1 + cos(MY_PI * dotpath));
    
    prec_t prefactor = 0.0;

    if (system_type == CLIMBER)
        prefactor = -2.0 * sys_info[sys_id].dot;
    else {
        prefactor = -sys_info[sys_id].dot + k_spring * (sys_info[sys_id].nlen - sys_info[sys_id].plen);
    }

    for (int i = 0; i < atom_num; i++)
    {
        f[i][0] += prefactor * tangent[i][0] +
            AngularContr * (springF[i][0] - sys_info[sys_id].dotSpringTangent * tangent[i][0]);
        f[i][1] += prefactor * tangent[i][1] +
            AngularContr * (springF[i][1] - sys_info[sys_id].dotSpringTangent * tangent[i][1]);
        f[i][2] += prefactor * tangent[i][2] +
            AngularContr * (springF[i][2] - sys_info[sys_id].dotSpringTangent * tangent[i][2]);
    }
}

void KMC::neb::out_stru(std::string temp_pos_path) const
//void KMC::neb::out_stru(const string &header) const
{
    //minimize->io_p->print_stru(*initial_system, header+to_string(0));
    for (int i = 0 ; i < n_image - 2; i++)
    {
        //string pos_name = header + to_string(i+1);
        //minimize->io_p->print_stru(inner_system[i], pos_name);
        string pos_name = temp_pos_path + "/POSCAR_" + to_string(i+1);
	minimize->io_p->print_stru(inner_system[i], pos_name);
    }
    //minimize->io_p->print_stru(*final_system, header+to_string(n_image-1));
}
