#include "potential.h"
using namespace KMC;
using namespace std;

KMC::ml_potential::ml_potential():potential(){}

KMC::ml_potential::ml_potential(const std::string &model):potential()
{
    dp.init(model);
}

void KMC::ml_potential::init(const std::string &model)
{
    dp.init(model);
}

prec_t KMC::ml_potential::infer(cell *system)
{
    prec_t energy;
    vector<prec_t> forces;
    vector<prec_t> virials;
    vector<prec_t> box(9);
    for (int i = 0 ; i < 9 ; ++i)
        box[i] = system->get_box()[i];
    dp.compute(energy, forces, virials, system->positions.get_coords(), system->get_atype(), box);
    system->forces.set_coords(forces);
    return energy;
}

vector<prec_t> KMC::ml_potential::infer(vector<cell> &systems)
{
    vector<prec_t> energies;
    vector<prec_t> forces;
    vector<prec_t> virials;
    vector<prec_t> boxs;
    vector<prec_t> positions;
    
    vector<prec_t> box(9);
    for (int i = 0 ; i < 9 ; ++i)
        box[i] = systems[0].get_box()[i];
    
    for (auto &sys : systems)
    {
        boxs.insert(boxs.end(), box.begin(), box.end());
        const vector<prec_t> temp_position = sys.positions.get_coords();
        positions.insert(positions.end(), temp_position.begin(), temp_position.end());
    }
    dp.compute(energies, forces, virials, positions, systems[0].get_atype(), boxs);

    size_t atom_num = systems[0].get_atom_num();
    vector<prec_t> system_forces;
    for (int i = 0 ; i < systems.size() ; ++i)
    {
        system_forces.assign(forces.begin()+i*atom_num*3, forces.begin()+(i+1)*atom_num*3);
        systems[i].forces.set_coords(system_forces);
    }
    return energies;
}