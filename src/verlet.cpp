#include "verlet.h"
#include <random>
#include <fstream>
#include "deepmd/DeepPot.h"
#include <iomanip> 

void KMC::Verlet::conf_md(string root_dir,string subiniconf_dir, int num_step, int dump_freq, double temp)
{
    string ini_file = subiniconf_dir + "/ini_conf";
    opt_cell = io_obj.read_stru(ini_file);
    string head = io_obj.header;
    prec_t scale_factor = io_obj.scaling;
    vector<string> elem = io_obj.elements;
    int natoms = opt_cell.get_atom_num();
    const vector<prec_t> &mass = opt_cell.get_mass();
    const atom_vec<bool_t> &moves = opt_cell.get_moves();
    array<double,9> box = opt_cell.get_box();
    vector<index_t> atype = opt_cell.get_atype();
    vector<prec_t> coords = opt_cell.positions.get_coords();
    vector<double> ini_coord,new_coord;
    vector<double> x_old,y_old,z_old,x_new,y_new,z_new,x_ini,y_ini,z_ini;
    vector<double> fx_old,fy_old,fz_old,fx_new,fy_new,fz_new;
    vector<double> vx_old,vy_old,vz_old,vx_new,vy_new,vz_new;
    x_ini.resize(natoms), y_ini.resize(natoms), z_ini.resize(natoms);
    x_old.resize(natoms), y_old.resize(natoms), z_old.resize(natoms);
    x_new.resize(natoms), y_new.resize(natoms), z_new.resize(natoms);
    fx_old.resize(natoms), fy_old.resize(natoms), fz_old.resize(natoms);
    fx_new.resize(natoms), fy_new.resize(natoms), fz_new.resize(natoms);
    vx_old.resize(natoms), vy_old.resize(natoms), vz_old.resize(natoms);
    vx_new.resize(natoms), vy_new.resize(natoms), vz_new.resize(natoms);

    string pot_path = root_dir + "/graph.000.pb";
    deepmd::DeepPot dp(pot_path);
    
    vector<double> tmp_box;
    tmp_box.clear();
    for (int i = 0 ; i < box.size(); i++)
       tmp_box.push_back(box[i]);     
    double lx = box[0], ly = box[4], lz = box[8];	       
    // initialize velocities randomly    
    vector<double> atomic_mass(natoms);
    vector<vector<double>> velocities(natoms);
    int n_atom = 0;
    unordered_map<index_t, size_t> ele_num;
    for (int i = 0 ; i < elem.size() ; i++)
        ele_num.insert(pair<index_t, size_t>(i, 0));
    for (auto type : atype)
        ele_num[type] += 1;
    for (int i = 0 ; i < elem.size(); i++)
    {	
        for (int j = 0 ; j < ele_num[i]; j++)
        {
            n_atom += 1;
            atomic_mass[n_atom-1] = mass[i];
        }
    }
    temp = 300.0;
    for (int i = 0; i < natoms; i++) 
    {   
        velocities[i] = gen_velocity(temp,atomic_mass[i]);
        vx_old[i] = velocities[i][0];
        vy_old[i] = velocities[i][1];
        vz_old[i] = velocities[i][2];

        x_old[i] = coords[3*i];
        y_old[i] = coords[3*i+1];
        z_old[i] = coords[3*i+2];
    }       
    auto maxElement = max_element(z_old.begin(), z_old.end());
    double max_z = *maxElement + 1.0;
    //The force at t = 0
    string log_path = subiniconf_dir + "/energy_md.log";
    ofstream ofs(log_path,ios::out);
    double e;
    ini_coord = coords;
    vector <double> f_old(3*natoms),f_new(3*natoms),v(3*natoms);
    dp.compute (e, f_old, v, ini_coord, atype, tmp_box);
    ofs << setprecision(9) << setw(16) << "0" << setw(16) << e << endl;
    for (int i = 0; i < natoms; i++)
    {
        fx_old[i] = f_old[3*i];
        fy_old[i] = f_old[3*i+1];
        fz_old[i] = f_old[3*i+2];
    }
    //dump the coordinate
    for (int step = 1; step <= num_step; step++)
    {
        if (step % dump_freq == 0)
        {
            string tmp_path = subiniconf_dir + "/POSCAR_" + to_string(step);
            io_obj.print_stru(opt_cell, tmp_path);
        }
        //update the positions using the velocity-Verlet algorithm  
        for (int i = 0; i < natoms; i++)
        {
            if (moves[i][0] == 1 && moves[i][1] == 1 && moves[i][2] == 1)
            { 
                x_new[i] = x_old[i] + vx_old[i]*dt + 0.5*fx_old[i]*dt*dt / atomic_mass[i]*convert_factor;
                y_new[i] = y_old[i] + vy_old[i]*dt + 0.5*fy_old[i]*dt*dt / atomic_mass[i]*convert_factor;
                z_new[i] = z_old[i] + vz_old[i]*dt + 0.5*fz_old[i]*dt*dt / atomic_mass[i]*convert_factor;
                if (x_new[i] < 0.0) x_new[i] = x_new[i] + lx;
                if (y_new[i] < 0.0) y_new[i] = y_new[i] + ly;
                if (z_new[i] < 0.0) z_new[i] = z_new[i] + lz;
                if (x_new[i] > lx) x_new[i] = x_new[i] - lx;
                if (y_new[i] > ly) y_new[i] = y_new[i] - ly;
                if (z_new[i] > lz) z_new[i] = z_new[i] - lz;
            }
            else
            {
                x_new[i] = x_old[i];
                y_new[i] = y_old[i];
                z_new[i] = z_old[i];
            }
        }
        // compute force at various timestpes
	new_coord.clear();
	for (int i = 0; i < natoms; i++)
        {
	    new_coord.push_back(x_new[i]);
            new_coord.push_back(y_new[i]);
            new_coord.push_back(z_new[i]);
	}
        dp.compute (e, f_new, v, new_coord, atype, tmp_box);
	ofs << setprecision(9) << setw(16) << step << setw(16) << e << endl;
	vector<double> force_wall; 
	force_wall = lj_force(z_new,max_z,natoms);
        for (int i = 0; i < natoms; i++)
        {
            int num_x = 3*i;
            int num_y = 3*i + 1;
            int num_z = 3*i + 2;
            fx_new[i] = f_new[num_x];
            fy_new[i] = f_new[num_y];
            fz_new[i] = f_new[num_z]+force_wall[i];
        }
        // update velocities using Verlet algorithm
        for (int i = 0; i < natoms; i++) 
        {
            if (moves[i][0] == 1 && moves[i][1] == 1 && moves[i][2] == 1)  
            {  
               vx_new[i] = vx_old[i] + 0.5*(fx_old[i] + fx_new[i])*dt / atomic_mass[i]*convert_factor;
               vy_new[i] = vy_old[i] + 0.5*(fy_old[i] + fy_new[i])*dt / atomic_mass[i]*convert_factor;
               vz_new[i] = vz_old[i] + 0.5*(fz_old[i] + fz_new[i])*dt / atomic_mass[i]*convert_factor;
            }
            else
            {
               vx_new[i] = 0.0;
               vy_new[i] = 0.0;
               vz_new[i] = 0.0;
            }        
            fx_old[i] = fx_new[i];
            fy_old[i] = fy_new[i];
            fz_old[i] = fz_new[i];

            vx_old[i] = vx_new[i];
            vy_old[i] = vy_new[i];
            vz_old[i] = vz_new[i];

            x_old[i] = x_new[i];
            y_old[i] = y_new[i];
            z_old[i] = z_new[i];
        }
        vector<prec_t> tmp_coords;
        tmp_coords.clear();
	for (int i = 0; i < natoms; i++)
        {
	    tmp_coords.push_back(x_new[i]);
	    tmp_coords.push_back(y_new[i]);
	    tmp_coords.push_back(z_new[i]);
	}
        opt_cell.positions.coords = tmp_coords;
	if (step % dump_freq == 0)
        {
            string tmp_path = subiniconf_dir + "/POSCAR_" + to_string(step);
            io_obj.print_stru(opt_cell, tmp_path);
        }
    }
}

vector<double> KMC::Verlet::gen_velocity(double temp, double mass) 
{
    random_device rd;
    mt19937 gen(rd());

    double standard_dev = sqrt(kB*temp/(mass/NA/1000));
    normal_distribution<double> dist(0.0, standard_dev);

    vector<double> velocities(3);
    velocities[0] = dist(gen)*1.0e-2; 
    velocities[1] = dist(gen)*1.0e-2; 
    velocities[2] = dist(gen)*1.0e-2; 

    return velocities;
}

vector<double> KMC::Verlet::lj_force(vector<double> z_coord,double z_wall, int natoms)
{
    double epsilon = 0.1;
    double sigma = 1.0;
    double cutoff = 2.0;
    vector<double> lj_force;
    lj_force.resize(natoms);
    for (int i = 0; i < natoms; i++)
    {
        double dz = z_coord[i] - z_wall;
	dz = sqrt(dz*dz);
	if (dz < cutoff)
        {
	    double distance = sigma / dz;
            lj_force[i] = -24*epsilon*(2*pow(distance, 12) - pow(distance, 6)) / dz;
	}
	else
	    lj_force[i] = 0.0;	
    }	    
    return lj_force;
}
