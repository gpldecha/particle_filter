#include <particle_filter/sampling.h>

namespace pf{


Sampling::Sampling(std::size_t D){

    A.resize(D,D);
    A.eye();
    z.resize(D);
    mu.resize(D);

    bFirst = true;

}

void Sampling::set_covariance(const arma::mat& covariance){
    A = chol(covariance);
}

void Sampling::resample(arma::mat& particles, arma::colvec &weights){

    resample_low_weights(particles,weights);
    //stratified(particles,weights);

}

void Sampling::resample_low_weights(arma::mat& particles, arma::colvec& weights){

    //std::cout<<   "Resample" << std::endl;

    int index       = U_index(generator);
    int nbParticles = particles.n_rows;
    double max_w    =  arma::max(weights); // TODO no need to recompute this

    U_beta          = boost::uniform_real<float>(0,2*max_w);
    arma::mat X_tmp = particles;



    float beta = 0;
    for(int i = 0; i < nbParticles;i++){
       if(weights(i) < 0.9*max_w){
//        if(weights(i) == 0){

            beta = beta + U_beta(fib_generator);

            while(beta > weights(index)){
                beta = beta - weights(index);
                index = (index + 1) % nbParticles;
            }

            //std::cout<< "index: " << index << " i " <<i << " max_i: " << particles.n_rows <<std::endl;


            z.randn();
            particles.row(i) = X_tmp.row(index) +  (A * z).st();

       }
    }

    for(int i = 0; i < nbParticles;i++){
        weights(i) = 1/(double)nbParticles;
    }
   // std::cout<< "finished" << std::endl;


}


void Sampling::stratified(arma::mat& particles, arma::colvec& weights){

    arma::vec Q = arma::cumsum(weights);
    particles_tmp = particles;//resize(particles.n_rows,particles.n_cols);

    T = arma::randu<arma::vec>(weights.n_elem)/((double)weights.n_elem);


    for(std::size_t i = 0; i < T.n_elem;i++){
        T(i) = T(i) + ((double)i)/((double)weights.n_elem);
    }
    T(T.n_elem-1) = 1;

    std::size_t i =0;
    std::size_t j =0;


    index.clear();

    while(i < weights.n_elem){

        if(j >= Q.n_rows){
            index.push_back(j-1);
            i = i + 1;
        }else{

            if(T(i) < Q(j)){
                index.push_back(j);
                i = i + 1;
            }else{
                j = j + 1;
            }
        }

    }

 //   std::cout<< "index.size(): " << index.size() << std::endl;

  //  index.print("index") << std::endl;
    assert(index.size() == particles_tmp.n_rows);
  //  std::cout<< "--2--" << std::endl;
    for(std::size_t i = 0; i < index.size();i++){
        z.randn();

        particles.row(i) = particles_tmp.row(index[i]) + (A * z).st();

    }

}

void Sampling::gaussian_kernel(arma::colvec &x){
    z.randn();
    x = x + A * z;
}

}
