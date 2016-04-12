#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <particle_filter/base_particle_filter.h>
#include <particle_filter/sampling.h>
#include <particle_filter/reguliser.h>
#include <ros/ros.h>
#include <array>

namespace pf {


class Particle_filter : public Base_particle_filter{

public:

    Particle_filter(const pf::likelihood_model& likelihood_function,
                    const pf::Measurement_h& measurement_h,
                    const pf::motion_model& motion_function,
                    std::size_t number_particles,
                    std::size_t x_dimension,
                    std::size_t y_dimension);

    virtual void update(const arma::colvec& u, const arma::colvec& Y, double duration=0.0)                   = 0;

    virtual void motion_update(const arma::colvec &u, double duration=0.0)             = 0;

    virtual void measurement_update(const arma::colvec& Y)                              = 0;

    virtual void init_visualise(ros::NodeHandle& node,const std::string& topic_name);

    virtual void compute_color(color_type c_type);

    virtual void visualise();

    void reinitialise(const arma::mat& points);

    virtual void set_rotation(const arma::mat33& Rot);

    void normalise();

    void reset_weights();

    void print();

    void add_reguliser(std::shared_ptr<pf::Reguliser>& reguliser);

public:

    arma::mat                    particles;
    arma::colvec                 weights;
    arma::colvec                 weights_tmp;
    arma::colvec                 L;
    arma::mat                    hY;
    double                       max_w;
    double                       sum_ww;
    double                       sample_variance;

protected:

    const pf::likelihood_model& likelihood_function;
    const pf::Measurement_h&    measurement_h;
    const pf::motion_model&     motion_function;
    std::size_t                 number_particles;
    double                      sum_w;
    double                       Y_dim;
    arma::mat33                 Rot;

    std::shared_ptr<pf::Reguliser>  reguliser;



};



class Particle_filter_sir : public Particle_filter {

public:

    Particle_filter_sir(const pf::likelihood_model& likelihood_function,
                        const pf::Measurement_h& measurement_h,
                        const pf::motion_model& motion_function,
                        std::size_t number_particles,
                        std::size_t x_dimension,
                        std::size_t y_dimension,
                        pf::Sampling& sampling);

    virtual void update(const arma::colvec &u, const arma::colvec &Y, double duration=0);

    virtual void motion_update(const arma::colvec &u, double duration=0);

    virtual void measurement_update(const arma::colvec& Y);


private:

    pf::Sampling&                   sampling;



};

}


#endif
