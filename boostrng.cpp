#include "boostrng.h"


#ifdef _BOOSTRNG_CPP_
#undef _BOSSTRNG_CPP_

BoostRNG::BoostRNG(uint _seed) : seed(_seed), rng(boost::mt19937(_seed)),
                                 uniformDistribution(boost::random::uniform_real_distribution<> (1e-6, 1.)),
                                 generateRandomNumbers(boost::variate_generator<boost::mt19937&,boost::random::uniform_real_distribution<> > (rng, uniformDistribution)){
 for(int i = 0; i < 1e6; i++) generateRandomNumbers();
}

double BoostRNG::get_number() {
    return generateRandomNumbers();
}

#endif

