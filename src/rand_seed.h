#ifndef RAND_H
#define RAND_H

#include <random>

using namespace std;

class RandSeed 
{
    public:
        RandSeed();
    public:
        double get_rand_seed();
    public:
        random_device rd_;
        mt19937 gen_;
        uniform_real_distribution<> dis_;
};

#endif
