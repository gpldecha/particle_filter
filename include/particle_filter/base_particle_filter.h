#ifndef BASE_PARTICLE_FILTER_H_
#define BASE_PARTICLE_FILTER_H_

#include <armadillo>
#include <ros/ros.h>
#include <boost/scoped_ptr.hpp>
#include <array>
#include <visualise/vis_point_cloud.h>
#include <visualise/colormap.h>
#include <particle_filter/particle_filter_definitions.h>


namespace pf {

class Base_particle_filter{


public:

    virtual void update(const arma::colvec& u, const arma::colvec& Y, double duration=0)    = 0;

    virtual void motion_update(const arma::colvec &u,double duration=0.0)                   = 0;

    virtual void measurement_update(const arma::colvec& Y)                                  = 0;

    virtual void init_visualise(ros::NodeHandle& node,const std::string& topic_name)        = 0;

    virtual void compute_color(pf::color_type c_type=pf::C_WEIGHTS)                         = 0;

    virtual void visualise()                                                                = 0;

    void set_visualisation_mode(opti_rviz::Vis_point_cloud::display_mode mode);

    void set_color_mode(pf::color_type color_t);

    virtual void set_rotation(const arma::mat33& Rot)                                       = 0;


protected:

    boost::scoped_ptr<opti_rviz::Vis_point_cloud>   vis_pf;
    pf::color_type                                  color_t;


    std::vector<std::array<float,3> >       colors;
    unsigned char                           rgb[3];


};

}

#endif
