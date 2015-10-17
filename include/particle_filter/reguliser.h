#ifndef PF_REGULISER_H_
#define PF_REGULISER_H_

#include <armadillo>


namespace pf{

class Reguliser{

public:

    Reguliser(arma::colvec3& position, arma::mat33& rot);

    void update_position(const arma::colvec3& position,const arma::mat33& rot);
    ///
    /// \brief update, only called when resample is called
    /// \param Y -> [0-1,0-1], measurement
    /// \param x, eof position
    /// \param rot, eof orientation
    ///
    bool update(const arma::colvec& Y,arma::mat& particles,arma::colvec& weights);


public:

    arma::mat        covariance;

private:

    bool              bFirst;
    arma::colvec3     position;
    arma::mat33       rot;


};

}



#endif
