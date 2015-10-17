#include <iostream>
#include <particle_filter/particle_filter.h>
#include <particle_filter/particle_filter_gmm.h>
#include <functional>
#include <armadillo>



void likeli(arma::colvec& L, const arma::colvec& Y,const arma::mat& X){
    for(std::size_t i = 0; i < Y.n_elem;i++){
        L(i) = 1;
    }
}

void motion(arma::mat& X,const arma::colvec& u){

}


int main(int argc, char** argv)
{

    std::cout<< "== start test ==" << std::endl;


    pf::likelihood_model likelihood_f = std::bind(likeli,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    pf::motion_model     motion_f     = std::bind(motion,std::placeholders::_1,std::placeholders::_2);

    mean_shift::MeanShift_Parameters mean_shift_parameters;

    std::unique_ptr<pf::Particle_filter_gmm> particle_filter = std::unique_ptr<pf::Particle_filter_gmm>
            (
                    new pf::Particle_filter_gmm(likelihood_f,motion_f,100,3,mean_shift_parameters,20)
            );

    particle_filter->particles = arma::randu(100,3);
    particle_filter->particles.print("particles");


    //particle_filter->initialise_mean_shift(0.4);


    arma::colvec3 u;
    arma::colvec3 Y;

    std::cout<< "   motion update" << std::endl;
    particle_filter->motion_update(u);
    std::cout<< "   measurement update" << std::endl;
    particle_filter->measurement_update(Y);


    return 0;
}
