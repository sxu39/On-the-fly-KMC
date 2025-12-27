#ifndef KMC_COORD_H
#define KMC_COORD_H

#include <vector>
#include "type.h"
#include <exception>

namespace KMC
{
    template <class T>
    class atom_vec
    {
        public:
            std::vector<T> coords;
            size_t atom_num;
        public:
            atom_vec();
            explicit atom_vec(size_t);
            explicit atom_vec(const std::vector<T> &);
            atom_vec(const atom_vec &);
            virtual ~atom_vec() = default;
            void init_coords(size_t);
            void set_coords(const std::vector<T> &);
            void set_coords(T);
            inline size_t get_atom_num() const;
            inline const std::vector<T> & get_coords() const;
            operator bool() const;
            T *operator[](index_t);
            const T *operator[](index_t) const;
    };

    enum norm_style {MAX, INF, TWO};
    class coord : public atom_vec<prec_t>
    {
        public:
            coord() = default;
            explicit coord(size_t);
            explicit coord(const std::vector<prec_t> &);
            coord(const coord &);
            ~coord() = default;
            void init_coords(size_t);
            void set_coords(const std::vector<prec_t> &);
            void set_coords(prec_t value=0.0);
            prec_t norm_inf(const atom_vec<bool_t> &moves=atom_vec<bool_t>()) const;
            prec_t norm_max(const atom_vec<bool_t> &moves=atom_vec<bool_t>()) const;
            prec_t norm_sqr(const atom_vec<bool_t> &moves=atom_vec<bool_t>()) const;
    };
}

template <class T>
size_t KMC::atom_vec<T>::get_atom_num() const
{
    return atom_num;
}

template <class T>
const std::vector<T> & KMC::atom_vec<T>::get_coords() const
{
    return coords;
}

#endif
