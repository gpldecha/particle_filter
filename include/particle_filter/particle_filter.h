#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <particle_filter/particle_filter_definitions.h>
#include <particle_filter/sampling.h>
#include <particle_filter/reguliser.h>
#include <visualise/vis_point_cloud.h>
#include <ros/ros.h>
#include <array>

namespace pf {


class Particle_filter{

public:

    Particle_filter(const pf::likelihood_model& likelihood_function,
                    const pf::motion_model& motion_function,
                    std::size_t number_particles,
                    std::size_t dimension);

    void reinitialise(const arma::mat& points);

    virtual void update(const arma::colvec& u, const arma::colvec& Y) = 0;

    virtual void motion_update(const arma::colvec &u) = 0;

    virtual void measurement_update(const arma::colvec& Y) = 0;

    void set_rotation(const arma::mat33& Rot);

    void normalise();

    void reset_weights();

    virtual void init_visualise(ros::NodeHandle& node,const std::string& topic_name);

    virtual void visualise();

    void  set_visualisation_mode(opti_rviz::Vis_point_cloud::display_mode mode);

    void set_color_mode(pf::color_type color_t);

    void print();

    void add_reguliser(std::shared_ptr<pf::Reguliser>& reguliser);


private:

    void compute_color(pf::color_type c_type=pf::C_WEIGHTS);

public:

    arma::mat                    particles;
    arma::colvec                 weights;
    arma::colvec                 weights_tmp;
    arma::colvec                 L;
    double                       max_w;
    double                       sum_ww;
    double                       sample_variance;
    pf::color_type               color_t;

protected:

    const pf::likelihood_model& likelihood_function;
    const pf::motion_model&     motion_function;
    std::size_t                 number_particles;
    double                      sum_w;
    arma::mat33                 Rot;

    std::shared_ptr<pf::Reguliser>  reguliser;

    std::unique_ptr<opti_rviz::Vis_point_cloud> vis_pf;


private:

    std::vector<std::array<float,3> >       colors;
    unsigned char                           rgb[3];


};



class Particle_filter_sir : public Particle_filter {

public:

    Particle_filter_sir(const pf::likelihood_model& likelihood_function,
                        const pf::motion_model& motion_function,
                        std::size_t number_particles,
                        std::size_t dimension,
                        pf::Sampling& sampling);

    virtual void update(const arma::colvec &u, const arma::colvec &Y);

    virtual void motion_update(const arma::colvec &u);

    virtual void measurement_update(const arma::colvec& Y);


private:

    pf::Sampling&                   sampling;



};

}


#endif
