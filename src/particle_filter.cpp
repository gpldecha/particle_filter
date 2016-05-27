#include <particle_filter/particle_filter.h>
#include <limits>
#include <visualise/colormap.h>

namespace pf{


Particle_filter::Particle_filter(const likelihood_model &likelihood_function,
                                 const Measurement_h &measurement_h,
                                 const motion_model &motion_function,
                                 std::size_t number_particles,
                                 std::size_t x_dimension,
                                 std::size_t y_dimension):
    likelihood_function(likelihood_function),
    measurement_h(measurement_h),
    motion_function(motion_function),
    number_particles(number_particles),
    Y_dim(y_dimension)
{

    particles.resize(number_particles,x_dimension);
    weights.resize(number_particles);
    weights_tmp.resize(number_particles);
    colors.resize(number_particles);
    L.resize(number_particles);

    assert(y_dimension != 0);
    assert(x_dimension != 0);

    hY.resize(number_particles,Y_dim);

   for(std::size_t i = 0; i < number_particles;i++){
        weights(i)      = 1;
        L(i)            = 1;
   }

    normalise();
    weights_tmp  = weights;

    color_t     = pf::C_WEIGHTS;

}

void Particle_filter::reinitialise(const arma::mat& points){

    number_particles = points.n_rows;
    particles        = points;
    weights.resize(number_particles);
    weights_tmp.resize(number_particles);
    colors.resize(number_particles);
    L.resize(number_particles);
    hY.resize(number_particles,Y_dim);

   for(std::size_t i = 0; i < number_particles;i++){
        weights(i)      = 1;
        L(i)            = 1;
   }

    normalise();
    weights_tmp  = weights;
}

void Particle_filter::reset_weights(){
    for(std::size_t i = 0; i < number_particles;i++){
         weights(i)      = 1;
         weights_tmp(i)  = 1;
         L(i)            = 1;
    }
    normalise();
    weights_tmp     = weights;
}

void Particle_filter::set_rotation(const arma::mat33& Rot){
    this->Rot = Rot;
}

void Particle_filter::normalise(){
    sum_w = arma::sum(weights) + std::numeric_limits<double>::min();
    if(sum_w == 0){
        sum_ww =0;
        max_w = 0;
        sample_variance = 0;
    }else{
        weights = weights / sum_w;
        sum_ww  = arma::dot(weights,weights);
        max_w  = arma::max(weights);
        sample_variance = 1.0/sum_ww;
    }
}

void Particle_filter::print(){
    std::cout<< "particles: " << particles.n_rows << " x " << particles.n_cols << std::endl;
    std::cout<< "weights:   " << weights.n_elem << " x 1 " << std::endl;
}

void Particle_filter::add_reguliser(std::shared_ptr<pf::Reguliser>& reguliser){
    this->reguliser = reguliser;
}

void Particle_filter::init_visualise(ros::NodeHandle& node,const std::string& topic_name){
    vis_pf.reset(new opti_rviz::Vis_point_cloud(node,topic_name));
    vis_pf->initialise("world",particles);
    vis_pf->set_display_type(opti_rviz::Vis_point_cloud::ONLY_HIGH_WEIGHTS);
}

void Particle_filter::compute_color(color_type c_type){

    if(c_type == C_LIKE){
        for(std::size_t i = 0 ; i < number_particles;i++){
            ColorMap::jetColorMap(rgb,L(i),0,arma::max(L));
            colors[i][0]    = ((float)rgb[0])/255;
            colors[i][1]    = ((float)rgb[1])/255;
            colors[i][2]    = ((float)rgb[2])/255;
        }
    }else{
        for(std::size_t i = 0 ; i < number_particles;i++){
            ColorMap::jetColorMap(rgb,weights(i),0,max_w);
            colors[i][0]    = ((float)rgb[0])/255;
            colors[i][1]    = ((float)rgb[1])/255;
            colors[i][2]    = ((float)rgb[2])/255;
        }
    }
}


void Particle_filter::visualise(){
    compute_color(color_t);
    vis_pf->update(particles,colors,weights,0.8 * max_w);
    vis_pf->publish();
}


///
/// \brief Particle_filter_sir::Particle_filter_sir
/// \param likelihood_function
/// \param motion_function
/// \param number_particles
/// \param dimension
/// \param sampling
///

Particle_filter_sir::Particle_filter_sir(const pf::likelihood_model& likelihood_function,
                                         const pf::Measurement_h& measurement_h,
                                         const pf::motion_model&     motion_function,
                                               std::size_t           number_particles,
                                               std::size_t           x_dimension,
                                               std::size_t           y_dimension,
                                               pf::Sampling&         sampling):
       Particle_filter(likelihood_function,measurement_h,motion_function,number_particles,x_dimension,y_dimension),
       sampling(sampling)
{

}



void Particle_filter_sir::update(const arma::colvec &u, const arma::colvec &Y, double duration){
    // std::cout<< "pf motion" << std::endl;
     motion_update(u,duration);
   //  std::cout<< "pf measurement" << std::endl;
     measurement_update(Y);
   //  std::cout<< "weights (" << weights.n_elem << " x 1 )" << std::endl;
   //  std::cout<< "particl (" << particles.n_rows << " x " << particles.n_cols << ")" << std::endl;



    if(reguliser != NULL){
        if(reguliser->update(Y,particles,weights)){
            normalise();
        }
    }

    if(sample_variance < 0.90 * static_cast<double>(number_particles)){
       //sampling.set_covariance(reguliser.covariance);
       sampling.resample(particles,weights);
       normalise();
    }
}

void Particle_filter_sir::motion_update(const arma::colvec& u, double duration){
    motion_function(particles,u);
}

void Particle_filter_sir::measurement_update(const arma::colvec& Y){

 //   std::cout<< "before measurement_h" << std::endl;
    //std::cout<< "hY:        (" << hY.n_rows        << " x " << hY.n_cols << ") " << std::endl;
    //std::cout<< "particles: (" << particles.n_rows << " x " << particles.n_cols << std::endl;

    measurement_h(hY,particles,Rot);

   // std::cout<< "before likelihood_function" << std::endl;

    //likelihood_function(L.memptr(),Y,hY);

  //  std::cout<< "after likelihood" << std::endl;

   //  check they are not all zero
    if(arma::sum(L) == 0){
        std::cout << "Particle_filter_sir: likelihood is zero" << std::endl;
        L = arma::ones(number_particles) * 1.0/(double) number_particles;
    }

   // std::cout<< "after likelihood_function" << std::endl;
   // std::cout<< "L       ( " << L.n_elem << " x 1)" << std::endl;
   // std::cout<<" weights (" << weights.n_elem << " x 1)" << std::endl;
    weights = L;// % weights_tmp;
    normalise();
    weights_tmp = weights;

}

/*void Particle_filter_sir::resample(){


}*/



}
