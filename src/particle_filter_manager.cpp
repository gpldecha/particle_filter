#include <particle_filter/particle_filter_manager.h>

namespace pf {

Particle_filter_manager::Particle_filter_manager(Particle_filter &particle_filter):
    particle_filter(particle_filter)
{

}

void Particle_filter_manager::update(const arma::fcolvec& Y,const arma::fcolvec& u){

    particle_filter.motion_update(u);
    particle_filter.measurement_update(Y);
    particle_filter.resample();

}



}
