#include <particle_filter/particle_filter_gmm.h>

namespace pf{

Particle_filter_gmm::Particle_filter_gmm(const likelihood_model&    likelihood_function,
                                         const Measurement_h&       measurement_h,
                                         const motion_model&        motion_function,
                                         const size_t               number_particles,
                                         const size_t               x_dimension,
                                         const size_t               y_dimension,
                                         const mean_shift::MeanShift_Parameters &mean_shift_parameters,
                                         const std::size_t num_intial_points_meanshift):
    Particle_filter(likelihood_function,measurement_h,motion_function,number_particles,x_dimension,y_dimension),
    mean_shift_parameters(mean_shift_parameters),
    num_intial_points_meanshift(num_intial_points_meanshift)
{

    bFirst = true;
    centroids.resize(num_intial_points_meanshift,x_dimension);


}

void Particle_filter_gmm::update(const arma::colvec &u, const arma::colvec &Y){
    motion_update(u);
    measurement_update(Y);

   /* if(reguliser != NULL){
        if(reguliser->update(Y,particles,weights)){
            normalise();
        }
    }*/

    if(bFirst){
       // std::cout<< "initialise_gmm -in-" << std::endl;
        initialise_gmm();
       // std::cout<< "initialise_gmm -out-" << std::endl;
        bFirst = false;
    }

    if(sample_variance < 0.9 * static_cast<double>(number_particles)){
        gmm.sample(particles);
        normalise();
        measurement_update(Y);
        estimate_gmm();
    }
}

void Particle_filter_gmm::motion_update(const arma::colvec &u){
    motion_function(particles,u);

    for(std::size_t k = 0; k < gmm.K;k++){
       //gmm.setMu(gmm.getMu(k) + u,k);
    }

}

void Particle_filter_gmm::measurement_update(const arma::colvec &Y){
    /// Update the weights
    //likelihood_function(L,Y,particles,Rot);


   /* std::cout<<std::endl;
    for(std::size_t i = 0; i < 5;i++){
       std::cout<< "L("<<i<<"): " << L(i) << std::endl;
    }*/


    weights = L;// (L + std::numeric_limits<double>::min());

    /*std::cout<<std::endl;
    for(std::size_t i = 0; i < 5;i++){
       std::cout<< "weights_1("<<i<<"): " << weights(i) << std::endl;
    }*/

    normalise();

    /*std::cout<<std::endl;
    for(std::size_t i = 0; i < 5;i++){
       std::cout<< "weights_2("<<i<<"): " << weights(i) << std::endl;
    }*/
}

void Particle_filter_gmm::initialise_gmm(){

  //  std::cout << "initialise_gmm" << std::endl;

    arma::uvec q1 = arma::find(L > 0.9*arma::max(L));
    particles_hw.resize(q1.n_elem,3);
    L_hw.resize(q1.n_elem);
    for(std::size_t i = 0; i < particles_hw.n_rows;i++){
        particles_hw.row(i) = particles.row(q1(i));
        L_hw(i)             = L(q1(i));
    }
  //  particles_hw.print("particles_hw");

    initialise.findInitialisation(particles_hw,centroids);

  //  centroids.print("initial centroids");

    mean_shift_ptr = std::unique_ptr<mean_shift::MeanShift>(new mean_shift::MeanShift(particles_hw,mean_shift_parameters));
    mean_shift_ptr->set_initial_center_guess(centroids);
    mean_shift_ptr->update();
    mean_shift_ptr->merge_centroids(mean_shift_parameters.min_merge_distance);

    modes          = mean_shift_ptr->modes;

 //   modes.print("initial modes");

    modes = modes.st();

    clustering.setInitCentroids(modes);

    particles_hw   = particles_hw.st();
    clustering.kmeans(particles_hw);
    particles_hw   = particles_hw.st();
    clustering.mixture_model(gmm,particles_hw,L_hw);
}

void Particle_filter_gmm::estimate_gmm(){

    arma::uvec q1 = arma::find(L > 0.95*arma::max(L));
    particles_hw.resize(q1.n_elem,3);
    L_hw.resize(q1.n_elem);
    for(std::size_t i = 0; i < particles_hw.n_rows;i++){
        particles_hw.row(i) = particles.row(q1(i));
        L_hw(i)             = L(q1(i));
    }

    initialise.findInitialisation(particles_hw,centroids);
    mean_shift_ptr->set_initial_center_guess(centroids);
    mean_shift_ptr->update();
    mean_shift_ptr->merge_centroids(mean_shift_parameters.min_merge_distance);
    modes          = mean_shift_ptr->modes;

    modes = modes.st();

    clustering.setInitCentroids(modes);

    particles_hw   = particles_hw.st();
    clustering.kmeans(particles_hw);
    particles_hw   = particles_hw.st();

    clustering.mixture_model(gmm,particles_hw,L_hw);
}

void Particle_filter_gmm::init_visualise(ros::NodeHandle& node,const std::string& topic_name){
    Particle_filter::init_visualise(node,topic_name);
    vis_gmm = std::unique_ptr<opti_rviz::Vis_gmm>
             (
                new opti_rviz::Vis_gmm(node,topic_name + "_gmm")
             );

    vis_gmm->initialise("world_frame",gmm.pi,gmm.Means,gmm.Covariances);
}

void Particle_filter_gmm::visualise(){
    Particle_filter::visualise();
    /*if(gmm.gmm.Means().size() > 0){
        gmm.gmm.Means()[0].print("Mean[0]");
    }*/
    vis_gmm->update(gmm.pi,gmm.Means,gmm.Covariances);
    vis_gmm->publish();
}




}
