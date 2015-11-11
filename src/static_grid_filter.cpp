#include <particle_filter/static_grid_filter.h>
#include <tf/LinearMath/Matrix3x3.h>

namespace pf {

Static_grid_filter::Static_grid_filter(const pf::likelihood_model& likelihood_function,
                                       const Measurement_h &measurement_h,
                   const pf::motion_model& motion_function,
                         std::size_t number_particles,
                         std::size_t x_dimension
                        , std::size_t y_dimension)
    :Particle_filter(likelihood_function,measurement_h,motion_function,number_particles,x_dimension,y_dimension)
{
    bFirst = true;
}


void Static_grid_filter::update(const arma::colvec& u, const arma::colvec& Y){
    motion_update(u);
    measurement_update(Y);
}

void Static_grid_filter::motion_update(const arma::colvec &u){
  //  motion_function(particles,u);
}

void Static_grid_filter::measurement_update(const arma::colvec& Y){
  //  likelihood_function(L,Y,particles,Rot);
    weights = L;// + 2*std::numeric_limits<double>::min();
    normalise();
    /*weights = weights % weights_tmp;
    normalise();
    weights_tmp = weights;*/
}

/*
void Static_grid_filter::init_visualise(ros::NodeHandle& node,const std::string& topic_name){
    vis_grid = std::unique_ptr<opti_rviz::Vis_gird>
             (
                new opti_rviz::Vis_gird(node,topic_name)
             );

    vis_grid->initialise("/world",particles);
}*/
/*
void Static_grid_filter::visualise(){
    if(bFirst){
        alphas.resize(number_particles);
        bFirst = false;
    }

    for(std::size_t i = 0; i < weights.n_elem;i++){
            alphas(i) =   weights(i) / max_w;
    }

    std::cout<< "alphas.size(): " << alphas.n_elem << std::endl;
    std::cout<< "weights.n_elem:" << weights.n_elem << std::endl;
    std::cout<< "alphas(0):     " << alphas(0) << std::endl;
    std::cout<< "weights(0):    " << weights(0) << std::endl;
    std::cout<< "max_w:         " << max_w << std::endl;


    vis_grid->update(particles,alphas);
    vis_grid->publish();
}*/



void Static_grid_filter::create_cube(arma::mat& points,float l, float w, float h, float bin_w,
                                         const arma::colvec3 position, const arma::vec3& rpy){

    std::cout<< "-1" << std::endl;

    std::size_t num_l = l/bin_w;
    std::size_t num_w = w/bin_w;
    std::size_t num_h = h/bin_w;

    if(num_h == 0){
        num_h = 1;
    }
    if(num_w == 0){
        num_w = 1;
    }
    if(num_l == 0){
        num_l = 1;
    }

    tf::Matrix3x3 R;
    tf::Vector3 p;
    std::cout<< "-2" << std::endl;
    std::cout<< "num_l: " << num_l << std::endl;
    std::cout<< "num_w: " << num_w << std::endl;
    std::cout<< "nu_h:  " << num_h << std::endl;

    //R.setRPY(rpy(0),rpy(1),rpy(2));
    R.setIdentity();

    tf::Vector3 t(-l/2.0 + position(0),-w/2.0 + position(1),-h/2.0 + position(2));
    std::cout<< "-3" << std::endl;

    std::cout<< "before resize" << std::endl;
    std::cout<< "total num points: " << num_l * num_w * num_h << std::endl;
    points.resize(num_l * num_w * num_h,3);

    std::cout<< "after resize" << std::endl;

    std::size_t index = 0;
    float step = bin_w/2;

    for(std::size_t i = 0; i < num_l; i++){
        for(std::size_t j = 0; j < num_w; j++){
            for(std::size_t k=0; k < num_h; k++){

                p.setX((float)i * bin_w + step);
                p.setY((float)j * bin_w + step);
                p.setZ((float)k * bin_w + step);
                p = R * p + t;
                points(index,0) = p.x();
                points(index,1) = p.y();
                points(index,2) = p.z();
                index++;
            }
        }
    }

    std::cout<< "finished creating cube" << std::endl;


}



}
