#ifndef STATIC_GRID_FILTER_H_
#define STATIC_GRID_FILTER_H_

#include <particle_filter/particle_filter.h>
#include <visualise/vis_grid.h>

namespace pf {

class Static_grid_filter : public pf::Particle_filter {

public:

    Static_grid_filter(const pf::likelihood_model& likelihood_function,
                       const pf::Measurement_h& measurement_h,
                       const pf::motion_model& motion_function,
                             std::size_t number_particles,
                             std::size_t x_dimension,
                             std::size_t y_dimension
                            );

    virtual void update(const arma::colvec& u, const arma::colvec& Y);

    virtual void motion_update(const arma::colvec &u);

    virtual void measurement_update(const arma::colvec& Y);

    static void create_gaussian_cube(arma::mat& points,arma::colvec3& mu,float var, float l, float w);

    static void create_cube(arma::mat& points, float l, float w, float h, float bin_w,
                            const arma::colvec3 position, const arma::vec3 &rpy);

 /*   void init_visualise(ros::NodeHandle& node,const std::string& topic_name);

    void visualise();*/

private:


private:

     std::unique_ptr<opti_rviz::Vis_gird> vis_grid;
     arma::colvec   alphas;
     bool           bFirst;

};


}

#endif
