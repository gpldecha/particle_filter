#ifndef PARTICLE_FILTER_DEFINITIONS_H_
#define PARTICLE_FILTER_DEFINITIONS_H_

#include <functional>
#include <armadillo>

namespace pf{


typedef enum {C_LIKE,C_WEIGHTS} color_type;

typedef  std::function<void (arma::colvec& L, const arma::colvec& Y,const arma::mat& X, const arma::mat33& Rot)> likelihood_model;
typedef  std::function<void (arma::mat& X,const arma::colvec& u)> motion_model;

}

#endif
