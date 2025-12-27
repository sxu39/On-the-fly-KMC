#include "cell.h"
using namespace KMC;
using namespace std;

static constexpr double MAXIMGCOUNT = 16;

KMC::cell::cell(const array<prec_t, 9> &box, const vector<index_t> &atype, const vector<prec_t> &mass, const vector<prec_t> &positions, 
const vector<bool_t> &moves, pos_type p_type, move_type m_type):
box(box), atype(atype), mass(mass), atom_num(atype.size()), positions(positions), moves(moves), velocities(atom_num), forces(atom_num), 
p_type(p_type), m_type(m_type)
{
    if (atom_num != this->positions.get_atom_num()) throw "The atom number isn't consistent.";
}

KMC::cell::cell(const cell &Cell):box(Cell.box), atype(Cell.atype), mass(Cell.mass), atom_num(Cell.atom_num), positions(Cell.positions), moves(Cell.moves),
velocities(Cell.velocities), forces(Cell.forces), p_type(Cell.p_type), m_type(Cell.m_type){}

void KMC::cell::init(const cell &Cell)
{
    box = Cell.box;
    atype = Cell.atype;
    mass = Cell.mass;
    atom_num = Cell.atom_num;
    positions = Cell.positions;
    moves = Cell.moves;
    velocities = Cell.positions;
    forces = Cell.forces;
    p_type = Cell.p_type;
    m_type = Cell.m_type;
}

void KMC::cell::init(const array<prec_t, 9> &box, const vector<index_t> &atype, const vector<prec_t> &mass, const vector<prec_t> &positions, const vector<bool_t> &moves)
{
    this->box = box;
    this->atype = atype;
    this->mass = mass;
    atom_num = this->atype.size();
    this->positions.set_coords(positions);
    this->moves.set_coords(moves);
    this->velocities.init_coords(atom_num);
    this->forces.init_coords(atom_num);
    if (atom_num != this->positions.get_atom_num()) throw "The atom number isn't consistent.";
}

void KMC::cell::minimum_image(double &dx, double &dy, double &dz, bool triclinic) const
{
    if (triclinic == 0)
    {
        if (fabs(dx) > (MAXIMGCOUNT * box[0]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dx);
        while (fabs(dx) > 0.5 * box[0])
        {
            if (dx < 0.0) dx += box[0];
            else dx -= box[0];
        }
        if (fabs(dy) > (MAXIMGCOUNT * box[4]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dy);
        while (fabs(dy) > 0.5 * box[4])
        {
            if (dy < 0.0) dy += box[4];
            else dy -= box[4];
        }
        if (fabs(dz) > (MAXIMGCOUNT * box[8]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dz);
        while (fabs(dz) > 0.5 * box[8])
        {
            if (dz < 0.0) dz += box[8];
            else dz -= box[8];
        }
    } 
    else
    {
        if (fabs(dz) > (MAXIMGCOUNT * box[8]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dz);
        while (fabs(dz) > 0.5 * box[8])
        {
            if (dz < 0.0)
            {
                dz += box[8];
                dy += box[5]; // yz
                dx += box[2]; // xz
            }
            else
            {
                dz -= box[8];
                dy -= box[5]; // yz
                dx -= box[2]; // xz
            }
        }
        if (fabs(dy) > (MAXIMGCOUNT * box[4]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dy);
        while (fabs(dy) > 0.5 * box[4])
        {
            if (dy < 0.0) {
                dy += box[4];
                dx += box[1]; // xy
            }
            else
            {
                dy -= box[4];
                dx -= box[1]; // xy
            }
        }
        if (fabs(dx) > (MAXIMGCOUNT * box[0]))
            throw("Atoms have moved too far apart ({}) for minimum image\n", dx);
        while (fabs(dx) > 0.5 * box[0]) {
            if (dx < 0.0) dx += box[0];
            else dx -= box[0];
        }
    }
}
