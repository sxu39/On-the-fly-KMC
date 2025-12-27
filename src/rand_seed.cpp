#include "rand_seed.h"

RandSeed::RandSeed() : rd_(), gen_(rd_()), dis_(0.0, 1.0) {}

double RandSeed::get_rand_seed() 
{
    return dis_(gen_);
}
