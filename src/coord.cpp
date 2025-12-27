#include "coord.h"
using namespace KMC;
using namespace std;

template <class T>
KMC::atom_vec<T>::atom_vec():atom_num(0){}

template KMC::atom_vec<bool_t>::atom_vec();
template KMC::atom_vec<prec_t>::atom_vec();

template<class T>
KMC::atom_vec<T>::atom_vec(size_t atom_num):atom_num(atom_num)
{
    coords.resize(atom_num*3);
}

template <class T>
KMC::atom_vec<T>::atom_vec(const vector<T> &coords):coords(coords), atom_num(coords.size()/3)
{
    if (coords.size() % 3) throw "The number of atom coordinates is wrong.";
}

template KMC::atom_vec<bool_t>::atom_vec(const vector<bool_t> &);

template <class T>
KMC::atom_vec<T>::atom_vec(const atom_vec<T> &atomVec)
{
    coords = atomVec.coords;
    atom_num = atomVec.atom_num;
    if (atom_num * 3 != coords.size()) throw "The atom number is wrong.";
}

template KMC::atom_vec<bool_t>::atom_vec(const atom_vec<bool_t> &);

template <class T>
void KMC::atom_vec<T>::init_coords(size_t atom_num)
{
    coords.resize(atom_num*3);
}

template <class T>
void KMC::atom_vec<T>::set_coords(const vector<T> &coords)
{
    if (atom_num * 3 != coords.size()) throw "Invalid length of the vector.";
    for (int i = 0 ; i < atom_num * 3 ; ++i)
        this->coords[i] = coords[i];
}

template void KMC::atom_vec<bool_t>::set_coords(const vector<bool_t> &);

template <class T>
void KMC::atom_vec<T>::set_coords(T value)
{
    for (auto &crd : coords) crd = value;
}

template <class T>
KMC::atom_vec<T>::operator bool() const
{
    return atom_num ? true : false;
}

template <class T>
T *KMC::atom_vec<T>::operator[](index_t index)
{
    return &coords[3*index];
}

template prec_t *KMC::atom_vec<prec_t>::operator[](index_t);

template <class T>
const T *KMC::atom_vec<T>::operator[](index_t index) const
{
    return &coords[3*index];
}

template const bool_t *KMC::atom_vec<bool_t>::operator[](index_t) const;
template const prec_t *KMC::atom_vec<prec_t>::operator[](index_t) const;

KMC::coord::coord(size_t atom_num):atom_vec<prec_t>(atom_num)
{
    set_coords();
}

KMC::coord::coord(const vector<prec_t> &coords):atom_vec<prec_t>(coords){}

KMC::coord::coord(const coord &Coord):atom_vec<prec_t>(Coord){}

void KMC::coord::init_coords(size_t atom_num)
{
    atom_vec<prec_t>::init_coords(atom_num);
    set_coords();
}

void KMC::coord::set_coords(const vector<prec_t> &coords)
{
    atom_vec<prec_t>::set_coords(coords);    
}

void KMC::coord::set_coords(prec_t value)
{
    atom_vec<prec_t>::set_coords(value);
}

prec_t KMC::coord::norm_inf(const KMC::atom_vec<bool_t> &moves) const
{
    prec_t max_norm = 0;
    if (moves)
    {
        const vector<bool_t> &move = moves.get_coords();
        for (int i = 0 ; i < atom_num * 3 ; ++i)
            max_norm = max(max_norm, coords[i]*coords[i]*move[i]);
    }
    else
        for (auto &crd : coords)
            max_norm = max(max_norm, crd*crd);
    return max_norm;
}

prec_t KMC::coord::norm_max(const KMC::atom_vec<bool_t> &moves) const
{
    prec_t max_norm = 0;
    if (moves)
    {
        const vector<bool_t> &move = moves.get_coords();
        for (int i = 0 ; i < atom_num ; i+=3)
            max_norm = max(max_norm, coords[i]*coords[i]*move[i]+coords[i+1]*coords[i+1]*move[i+1]+coords[i+2]*coords[i+2]*move[i+2]);
    }
    else
        for (int i = 0 ; i < atom_num ; i+=3)
            max_norm = max(max_norm, coords[i]*coords[i]+coords[i+1]*coords[i+1]+coords[i+2]*coords[i+2]);
    return max_norm;
}

prec_t KMC::coord::norm_sqr(const KMC::atom_vec<bool_t> &moves) const
{
    prec_t norm = 0;
    if (moves)
    {
        const vector<bool_t> &move = moves.get_coords();
        for (int i = 0 ; i < atom_num * 3 ; ++i)
            norm += coords[i]*coords[i]*move[i];
    }
    else
        for (auto &crd : coords)
            norm += crd*crd;
    return norm;
}
