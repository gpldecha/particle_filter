#ifndef PARTICLE_FILTER_SAMPLING_H_
#define PARTICLE_FILTER_SAMPLING_H_

#include <armadillo>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace pf{

class Sampling_Parameters{
public:

    Sampling_Parameters(){
        kernel_covariance.zeros();
        kernel_covariance(0,0) = std::pow(0.0005,2);
        kernel_covariance(1,1) = std::pow(0.002,2);
        kernel_covariance(2,2) = std::pow(0.002,2);
    }

    void set_diag(double sd_x, double sd_y, double sd_z){
        kernel_covariance(0,0) = sd_x * sd_x;
        kernel_covariance(1,1) = sd_y * sd_y;
        kernel_covariance(2,2) = sd_z * sd_z;
    }

    arma::mat33 kernel_covariance;

};


class Sampling{

    typedef boost::mt19937                              ENG;              // Mersenne Twister
    typedef boost::normal_distribution<float>          GAUSSIAN_DIST;    // Normal Distribution
    typedef boost::uniform_real<float>                 UNIFORM_REAL;     // uniform real distribution
    typedef boost::uniform_int<int>                     UNIFORM_DISCRETE; // uniform discrete distribution
    typedef boost::variate_generator<ENG,GAUSSIAN_DIST> GAUSSIAN_GEN;


public:

    Sampling(std::size_t D);

    void resample(arma::mat &particles, arma::colvec &weights);

    void set_covariance(const arma::mat& covariance);

private:

    void resample_low_weights(arma::mat &particles, arma::colvec &weights);

    void stratified(arma::mat& particles, arma::colvec& weights);

    void gaussian_kernel(arma::colvec& x);

private:

    ENG                         generator;
    GAUSSIAN_DIST               gaussian_dist;
    UNIFORM_REAL                U_beta;
    UNIFORM_DISCRETE            U_index;
    boost::lagged_fibonacci607  fib_generator;
    bool                        bFirst;

    arma::colvec               noise;
    arma::mat                  covariance;
    arma::mat                  A;
    arma::colvec               z;
    arma::colvec               mu;
    arma::colvec                T;
    std::vector<unsigned int>   index;
    arma::mat                   particles_tmp;
};

}


#endif
