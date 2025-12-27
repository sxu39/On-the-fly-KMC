#ifndef KMC_CELL_H
#define KMC_CELL_H

#include <array>
#include <cmath>
#include "coord.h"

namespace KMC
{
    enum cell_type {SINGLE, INITIAL, INNER, FINAL, CLIMBER};
    enum pos_type {CARTESIAN, DIRECT};
    enum move_type {FREE, SELECTIVE};

    class cell
    {
        public:
            std::array<prec_t, 9> box;
            std::vector<index_t> atype;
            std::vector<prec_t> mass;
            size_t atom_num;
            atom_vec<bool_t> moves;
        public:
            coord positions;
            coord velocities;
            coord forces;
            pos_type p_type;
            move_type m_type;
        private:
            friend class neb;
            void minimum_image(double &, double &, double &, bool triclinic=true) const;
        public:
            cell() = default;
            cell(const std::array<prec_t, 9> &, const std::vector<index_t> &, const std::vector<prec_t> &, const std::vector<prec_t> &, 
                 const std::vector<bool_t> &, pos_type p_type=CARTESIAN, move_type m_type=FREE);
            cell(const cell &);
            ~cell() = default;
            void init (const cell &);
            void init(const std::array<prec_t, 9> &, const std::vector<index_t> &, const std::vector<prec_t> &, const std::vector<prec_t> &, const std::vector<bool_t> &);
            inline const std::array<prec_t, 9> & get_box() const;
            inline const std::vector<index_t> & get_atype() const;
            inline const std::vector<prec_t> & get_mass() const;
            inline const atom_vec<bool_t> & get_moves() const;
            inline size_t get_atom_num() const;
            friend class min;
    };
}

const std::array<KMC::prec_t, 9> & KMC::cell::get_box() const
{
    return box;
}

const std::vector<KMC::index_t> & KMC::cell::get_atype() const
{
    return atype;
}

const std::vector<KMC::prec_t> & KMC::cell::get_mass() const
{
    return mass;
}

const KMC::atom_vec<KMC::bool_t> & KMC::cell::get_moves() const
{
    return moves;
}

size_t KMC::cell::get_atom_num() const
{
    return atom_num;
}

#endif
