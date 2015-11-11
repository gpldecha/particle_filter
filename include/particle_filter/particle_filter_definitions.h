#ifndef PARTICLE_FILTER_DEFINITIONS_H_
#define PARTICLE_FILTER_DEFINITIONS_H_

#include <functional>
#include <armadillo>

namespace pf{


typedef enum {C_LIKE,C_WEIGHTS} color_type;

typedef  std::function<void (arma::colvec& L, const arma::colvec& Y, const arma::mat& hY)> likelihood_model;



typedef  std::function<void (arma::mat& X,const arma::colvec& u)>                             motion_model;

/**
 * @brief measurement_h : measurement function h(X) prototype, Y = h(X)
 *                        given state information X, produces a virtual sensation.
 */
typedef  std::function<void (arma::mat& hY,const arma::mat& points, const arma::mat33& Rot) > Measurement_h;


}

#endif
