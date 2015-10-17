#ifndef PARTICLE_FILTER_VISUALISATION_H_
#define PARTICLE_FILTER_VISUALISATION_H_

#include <ros/ros.h>

#include <armadillo>
#include <optitrack_rviz/include/visualise/vis_points.h>

namespace pf {

class ParticleFilterVis{

public:

    ParticleFilterVis(ros::NodeHandle& node, const std::string& topic_name):
        vis_points(node,topic_name)
    {

    }

    void initialise(const std::string& frame_id,const arma::mat& points){
        vis_points.initialise(frame_id,points);
    }

    void update(const arma::mat& points){
        vis_points.update(points);
    }

    void publish(){
        vis_points.publish();
    }


private:

    opti_rviz::Vis_points vis_points;

};


}

#endif
