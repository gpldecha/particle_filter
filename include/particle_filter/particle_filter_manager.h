#ifndef PARTICLE_FILTER_MANAGER_H_
#define PARTICLE_FILTER_MANAGER_H_

/**
    Particle filter manager base class
    ----------------------------------

    o Initialise Particle filter
    o Initialise prior probability distribution of particles
    o Initialise process and measurement noise

  **/

#include <distributions.h>
#include <particle_filter/particle_filter.h>


namespace pf{

class Particle_filter_manager{

public:

    Particle_filter_manager(Particle_filter&  particle_filter);

    virtual void initialise_prior_pdf() = 0;

    void update(const arma::fcolvec &Y, const arma::fcolvec& u);

protected:

    Distribution*      prior_pdf;
    Particle_filter&   particle_filter;

};


}


#endif
