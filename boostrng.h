# include <boost/random/mersenne_twister.hpp>
# include <boost/random/uniform_real.hpp>
# include <boost/random/variate_generator.hpp>

#ifndef _BOOSTRNG_H_
#define _BOOSTRNG_H_
#define _BOOSTRNG_CPP_

class BoostRNG {
  private:
    uint seed;
    boost::mt19937 rng;
    boost::random::uniform_real_distribution< > uniformDistribution;
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
        generateRandomNumbers;

  public:
    BoostRNG(uint seed);
    double get_number();
};

#endif
