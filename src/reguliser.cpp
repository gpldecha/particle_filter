#include <particle_filter/reguliser.h>
#include <statistics/distributions/distributions.h>

namespace pf {


Reguliser::Reguliser(arma::colvec3& position, arma::mat33& rot):
position(position),rot(rot){

    bFirst = true;
    covariance.resize(3,3);
    covariance.eye(3,3);
}

void Reguliser::update_position(const arma::colvec3& position,const arma::mat33& rot){
    this->position = position;
    this->rot = rot;
}

bool Reguliser::update(const arma::colvec& Y, arma::mat &particles, arma::colvec &weights){


    if(bFirst){
        // surface contact
        if(Y(0) == 1){
            std::cout<< " First Contact " << std::endl;
            // get the the boundaries of the weights which have positive weights
            double max_w = arma::max(weights);
            arma::uvec idx = arma::find(weights > 0.7*max_w);
            arma::mat X(idx.n_elem,3);
            for(std::size_t i = 0; i < X.n_rows;i++){
                X.row(i) = particles.row(idx(i));
            }

            std::cout<< "here" << std::endl;

            std::cout<< "X: " << X.n_rows << " x " << X.n_cols << std::endl;

            arma::mat max_elem = arma::max(X);
            arma::mat min_elem = arma::min(X);


            std::cout<< "here 2" << std::endl;

            max_elem.print("max_elem");
            min_elem.print("min_elem");

            arma::vec3 origin;
            origin(0) = (max_elem(0) - min_elem(0))/2;
            origin(1) = (max_elem(1) - min_elem(1))/2;
            origin(2) = (max_elem(2) - min_elem(2))/2;
            arma::mat orient;
            orient.eye(3,3);



            std::cout<< "here 3" << std::endl;

            double length = max_elem(0) - min_elem(0);
            double width  = max_elem(1) - min_elem(1);
            double height = max_elem(2) - min_elem(2);

            length = 0.01;
            width  = 0.8;
            height = 0.05;//height + 0.2*length;


            origin(0) = position(0);
            origin(1) = 0;
            origin(2) = position(2);


            origin.print("origin");
            std::cout<< "length: " << length << std::endl;
            std::cout<< "width: "  << width << std::endl;
            std::cout<< "height: " << height << std::endl;


            stats::Uniform uniform(origin,orient,length,width,height);
            std::cout<< "here" << std::endl;
            for(std::size_t i = 0; i < particles.n_rows;i++){
                particles.row(i) = uniform.sample().st();
            }

            weights.ones();
            weights = weights / (double) weights.n_elem;

            bFirst = false;
            return true;
       }else{
            covariance(0,0) = std::pow(0.02,2);
            covariance(1,1) = std::pow(0.02,2);
            covariance(2,2) = std::pow(0.02,2);
        }

    }else{

      /*  if(Y(0) == 1){

            covariance(0,0) = std::pow(0.0001,2);
            covariance(1,1) = std::pow(0.001,2);
            covariance(2,2) = std::pow(0.001,2);

        }else{

            covariance(0,0) = std::pow(0.02,2);
            covariance(1,1) = std::pow(0.02,2);
            covariance(2,2) = std::pow(0.02,2);


        }*/


    }

    return false;

}

}
