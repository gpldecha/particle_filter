#ifndef PARTICLE_FILTER_DEFINITIONS_H_
#define PARTICLE_FILTER_DEFINITIONS_H_

#include <functional>
#include <armadillo>
#include <tf/transform_broadcaster.h>

namespace pf{


typedef enum {C_LIKE,C_WEIGHTS} color_type;

//typedef  std::function<void (arma::colvec& L, const arma::colvec& Y, const arma::mat& hY)> likelihood_model;
typedef  std::function<void (double* L, const arma::colvec& Y, const arma::mat& hY)>       likelihood_model;



typedef  std::function<void (arma::mat& X,const arma::colvec& u)>                             motion_model;

/**
 * @brief measurement_h : measurement function h(X) prototype, Y = h(X)
 *                        given state information X, produces a virtual sensation.
 */
typedef  std::function<void (arma::mat& hY,const arma::mat& points, const arma::mat33& Rot) > Measurement_h;

//https://gerardus.googlecode.com/svn/trunk/matlab/GerardusCommon.hxx
template<typename T>
inline void ind2sub(T& i,
                    T& j,
                    T& k,
                    const int n_rows,
                    const int n_cols,
                    int n)
{

    // convert linear index to r, c, s
    k   =  n / (n_rows * n_cols); // slice value (Note: integer division)
    j   = (n % (n_rows * n_cols)) / n_rows; // column value (Note: integer division)
    i   = (n % (n_rows * n_cols)) % n_rows; // row value
}

inline std::size_t sub2ind_col_major(const std::size_t i,
                                     const std::size_t j,
                                     const std::size_t k,
                                     const std::size_t n_rows,
                                     const std::size_t n_cols){
    return i + n_rows * (j + n_cols * k);
}

template<typename T>
inline void tf_debuf(const arma::Col<T>& pos,const std::string& name){
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    transform.setOrigin( tf::Vector3(pos(0),pos(1),pos(2)) );
    tf::Quaternion q;
    q.setEuler(0,0,0);
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", name));
}


inline void tf_debuf(const arma::vec& pos,const std::string& name){
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    transform.setOrigin( tf::Vector3(pos(0),pos(1),pos(2)) );
    tf::Quaternion q;
    q.setEuler(0,0,0);
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", name));
}


inline void tf_debuf(const tf::Vector3& pos, const tf::Matrix3x3& Rot,const std::string& name){
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    transform.setOrigin( pos );
    tf::Quaternion q;
    Rot.getRotation(q);
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", name));
}





}

#endif
