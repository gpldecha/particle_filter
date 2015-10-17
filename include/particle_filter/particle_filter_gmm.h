#ifndef PARTICLE_FILTER_GMM_H_
#define PARTICLE_FILTER_GMM_H_

#include <particle_filter/particle_filter.h>

#include <statistics/meanshift.h>
#include <statistics/initialise.h>
#include <statistics/clustering.h>
#include <memory>
#include <visualise/vis_gmm.h>

namespace pf {



class Particle_filter_gmm : public Particle_filter{

public:

    Particle_filter_gmm(const pf::likelihood_model& likelihood_function,
                        const pf::motion_model& motion_function,
                        const std::size_t number_particles,
                        const std::size_t dimension,
                        const mean_shift::MeanShift_Parameters& mean_shift_parameters,
                        const std::size_t num_intial_points_meanshift);

    virtual void update(const arma::colvec &u, const arma::colvec &Y);

    virtual void motion_update(const arma::colvec &u);

    virtual void measurement_update(const arma::colvec& Y);

    virtual void init_visualise(ros::NodeHandle& node,const std::string& topic_name);

    virtual void visualise();

private:

    void initialise_gmm();

    void estimate_gmm();


public:

    GMM                                         gmm;

private:

    std::unique_ptr<mean_shift::MeanShift>      mean_shift_ptr;
    Clustering                                  clustering;
    Initialise                                  initialise;
    mean_shift::MeanShift_Parameters            mean_shift_parameters;
    std::size_t                                 num_intial_points_meanshift;

    arma::mat                                   centroids;
    arma::mat                                   modes;
    arma::mat                                   particles_hw;
    arma::vec                                   L_hw;
    std::vector<arma::vec>                      means;
    std::vector<arma::mat>                      covariances;

    std::unique_ptr<opti_rviz::Vis_gmm>         vis_gmm;
    bool                                        bFirst;
    double                                      numerical_precision;



};

}

#endif
