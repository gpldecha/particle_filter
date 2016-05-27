#ifndef POINT_MASS_FILTER_H_
#define POINT_MASS_FILTER_H_

#include "particle_filter/base_particle_filter.h"
#include "particle_filter/particle_filter_definitions.h"
#include <chrono>
#include <tf/transform_broadcaster.h>
#include <visualise/vis_gmm.h>

namespace pf {

class Point_mass_filter : public Base_particle_filter {

public:

    typedef struct delta{
    public:
        delta(float m=0.04, float n=0.04, float k=0.02, float min_delta=0.005):
            m(m),n(n),k(k),min_delta(min_delta){}
        float m,n,k,min_delta;
    } delta;

    typedef struct length{
    public:
        length(float m=0.5,float n=1.2,float k = 0.1):m(m),n(n),k(k){}
        float m,n,k;
    }length;


private:

    class BBox
    {
    public:

        BBox():max_x_n(0),
            min_x_n(0),
            max_y_n(0),
            min_y_n(0),
            max_z_n(0),
            min_z_n(0){}

        void print() const{
            std::cout<< "=== BBox ===" << std::endl;
            std::cout<< "x: " << max_x_n << "  " << min_x_n << std::endl;
            std::cout<< "y: " << max_y_n << "  " << min_y_n << std::endl;
            std::cout<< "z: " << max_z_n << "  " << min_z_n << std::endl;
        }

        void get_cropped_index(int& i_max_n,int& i_min_n,const int max_index,const int min_index,const float d,const int num_elem){
           // std::cout<< "get_cropped_index" << std::endl;
            int n = 0.02/d;
            if(n % 2 == 0){n++;}

            //std::cout<< "n: " << n << std::endl;

            int center_index = 0.5 * (max_index - min_index) + min_index;

            //std::cout<< "center_index: " << center_index << std::endl;

            i_min_n = center_index - static_cast<int>(0.5 * (n-1));
            i_max_n = center_index + static_cast<int>(0.5 * (n-1));

            if(i_max_n >= num_elem){
                i_max_n = num_elem;
            }
            if(i_min_n <= 0){
                i_min_n = 0;
            }

        }

        void update(const arma::mat& points, const arma::cube& P, const delta& delta_){
            tt(0) = points(max_x_n,0);
            tt(1) = points(max_y_n,1);
            tt(2) = points(max_z_n,2);

            bb(0) = points(min_x_n,0);
            bb(1) = points(min_y_n,1);
            bb(2) = points(min_z_n,2);

            dist_x = std::fabs(tt(0) - bb(0) + delta_.m);
            dist_y = std::fabs(tt(1) - bb(1) + delta_.n);
            dist_z = std::fabs(tt(2) - bb(2) + delta_.k);

            int i_max,j_max,kk_max,i_min,j_min,kk_min;

            pf::ind2sub<int>(i_max,j_max,kk_max,P.n_rows,P.n_cols,max_x_n);
            pf::ind2sub<int>(i_min,j_min,kk_min,P.n_rows,P.n_cols,min_x_n);
            if(dist_x < 0.01)
            {
                get_cropped_index(r_max,r_min,i_max,i_min,delta_.n,P.n_rows-1);
            }else{
                r_max = i_max;
                r_min = i_min;
            }

            pf::ind2sub<int>(i_max,j_max,kk_max,P.n_rows,P.n_cols,max_y_n);
            pf::ind2sub<int>(i_min,j_min,kk_min,P.n_rows,P.n_cols,min_y_n);
            if(dist_y < 0.05)
            {
                get_cropped_index(c_max,c_min,j_max,j_min,delta_.m,P.n_cols-1);
            }else{
                c_max = j_max;
                c_min = j_min;
            }

            pf::ind2sub<int>(i_max,j_max,kk_max,P.n_rows,P.n_cols,max_z_n);
            pf::ind2sub<int>(i_min,j_min,kk_min,P.n_rows,P.n_cols,min_z_n);
            if(dist_z < 0.05)
            {
                get_cropped_index(k_max,k_min,kk_max,kk_min,delta_.k,P.n_slices-1);
            }else{
                k_max = kk_max;
                k_min = kk_min;
            }


            int max_p_index = pf::sub2ind_col_major(r_max,c_max,k_max,P.n_rows,P.n_cols);
            int min_p_index = pf::sub2ind_col_major(r_min,c_min,k_min,P.n_rows,P.n_cols);

            tt = points.row(max_p_index).st();
            bb = points.row(min_p_index).st();

            dist_x = std::fabs(tt(0) - bb(0) + delta_.m);
            dist_y = std::fabs(tt(1) - bb(1) + delta_.n);
            dist_z = std::fabs(tt(2) - bb(2) + delta_.k);

        }

        /// indicies in column format of cube P
        std::size_t max_x_n, min_x_n;
        std::size_t max_y_n, min_y_n;
        std::size_t max_z_n, min_z_n;

        int r_max, r_min;
        int c_max, c_min;
        int k_max, k_min;

        /// indicies in (r,c,k) format of cube P

        arma::colvec3 tt,bb;
        double dist_x,dist_y,dist_z;


    };

    class GaussianKernel{
    public:

        GaussianKernel():var_x(1),var_y(1),var_z(1){}

        void reset(const delta& delta_){
            resize(var_x,var_y,var_z,delta_);
        }

        void print() const{
            std::cout<< "=== Gaussian Kernel ===" << std::endl;
            std::cout<< "var: " << var_x << " " << var_y << " " << var_z << std::endl;
            kernel_x.print("kernel_x");
            kernel_y.print("kernel_y");
            kernel_z.print("kernel_z");
        }

    private:

        void resize(double var_x, double var_y, double var_z, const delta& delta_){
            kernel_x.resize(get_kernel_size(var_x,delta_.m));
            kernel_y.resize(get_kernel_size(var_y,delta_.n));
            kernel_z.resize(get_kernel_size(var_z,delta_.k));

            set_kernel(kernel_x,var_x,delta_.m);
            set_kernel(kernel_y,var_y,delta_.n);
            set_kernel(kernel_z,var_z,delta_.k);

            if(!kernel_x.is_finite()){
               make_identity(kernel_x);
            }
            if(!kernel_y.is_finite()){
               make_identity(kernel_y);
            }

            if(!kernel_z.is_finite()){
               make_identity(kernel_z);
            }

        }

        void make_identity(arma::vec& kernel){
            kernel.zeros();
            int mid = (kernel.n_elem-1)/2;
            kernel(mid) = 1;
        }

        int get_kernel_size(double var,float bins){
            int dm = 2.0 * 3.0 * sqrt(var) / static_cast<double>(bins);
            if(dm <= 2)
            {
                return 3;
            }else if( dm >= 11){
                return 11;
            }else if( dm % 2 != 1)
            {
                return dm + 1;
            }else{
                return dm;
            }
        }

        void set_kernel(arma::vec& kernel,double var,float bins){

            int mid = (kernel.n_elem - 1)/2;
            double r;

            for(int i = 0; i < static_cast<int>(kernel.n_elem);i++){
                r = static_cast<double>(std::abs(i - mid)) * bins;
                kernel(i) = (1.0/sqrt(2.0 * M_PI * var)) * exp(-1.0/(2.0 * var) * r*r);
            }
            kernel = kernel / arma::sum(kernel);

        }

    public:
        arma::vec kernel_x, kernel_y,kernel_z;
        double var_x, var_y, var_z;

    };


public:

    Point_mass_filter(const pf::likelihood_model& likelihood_function,
                      const Measurement_h &measurement_h,
                      delta delta_ , length length_, std::size_t Y_dim);

    virtual void update(const arma::colvec& u, const arma::colvec& Y, double duration=0);

    void update_lik_debug(const arma::colvec &u, const arma::colvec &Y, double duration=0);

    void update_debug(const arma::colvec &u, const arma::colvec &Y, double duration=0);

    void update_real(const arma::colvec &u, const arma::colvec &Y, double duration=0);

    virtual void motion_update(const arma::colvec &u, double duration=0.0);

    virtual void measurement_update(const arma::colvec& Y);

    virtual void init_visualise(ros::NodeHandle& node,const std::string& topic_name);

    virtual void compute_color(pf::color_type c_type=pf::C_WEIGHTS);

    virtual void visualise();

    void reset(const arma::colvec3 &center, delta &delta_, length &length_);

    virtual void set_rotation(const arma::mat33& Rot);

    void print_prob() const ;

    const delta& get_delta() const;


private:

    /**
     * @brief truncate : sets probabilites to zero for gride cells
     * which are below < epsilon /  N.
     * The function then normalises again
     **/
    /*  void truncate();*/

    void normalise();

    void get_coordiantes();

    void increase_density();

    void decrease_density();

    void increase_lenght(int dm, int dn, int dk);

    void interpolate();

    double mean_weight() const;

    /**
     * @brief transform_bbox Given the top and bottom corner of a bounding box, removes all non
     *        zero points and creates a new cube.
     */
    void transform_bbox();

    /**
     * @brief get_bbox : finds the top and bottom corner grid cells of the point mass filter
     *                   which are non-zero.
     */
    void get_bbox();

    void get_index_boundary(int& max, int& min,int n, int num_cells);

    void get_boundary_values(double &max_r, double &max_c, double &max_k);

    arma::uword  index_cp;

private:


    inline double get_num_if_increase() const{
        return  (2.0*static_cast<double>(m_)-1.0) * (2.0*static_cast<double>(n_)-1.0) * (2.0*static_cast<double>(k_)-1.0);
    }

    inline double get_num_if_decrease() const{
        return  ((static_cast<double>(m_)-1.0)/2.0) * ((static_cast<double>(n_)-1.0)/2.0) * ((static_cast<double>(k_)-1.0)/2.0);
    }

    inline void get_num_non_zero(){
        normalise();
        max_w = P.max();
        eps_threashold = eps_trunc * max_w;
        assert(std::isfinite(eps_threashold));

        float sum_P = arma::sum(arma::vectorise(P));
        if((int)(10000 * sum_P)   != 10000){
            std::cout<< "sum_P: " << sum_P << std::endl;
            std::cout<< "sum_P: " << (int)(10000 * sum_P)  << std::endl;
        }
        assert( (int)(10000 * sum_P) == 10000);
        total_num_points    = P.n_elem;
        num_non_zero_points = total_num_points;

        for(std::size_t n = 0; n < P.n_elem;n++){
            if(P.at(n) < eps_threashold){
                num_non_zero_points = num_non_zero_points - 1;
                P.at(n) =0;
            }
        }
        normalise();
        max_w = P.max();
        eps_threashold = eps_trunc * max_w;
        assert(std::isfinite(eps_threashold));

    }


    inline void size_max_P(double &dx, double &dy, double dz){

        arma::colvec3 zero = points.row(0).st(); // P(0,0,0)
        arma::colvec3 px = points.row(pf::sub2ind_col_major(P.n_rows-1,0,0,P.n_rows,P.n_cols)).st(); // P(end,0,0)
        arma::colvec3 py = points.row(pf::sub2ind_col_major(0,P.n_cols-1,0,P.n_rows,P.n_cols)).st(); // P(0,end,0)
        arma::colvec3 pz = points.row(pf::sub2ind_col_major(0,0,P.n_slices-1,P.n_rows,P.n_cols)).st(); // P(0,0,end)


        dx = std::sqrt(arma::sum(arma::pow(px - zero,2)));
        dy = std::sqrt(arma::sum(arma::pow(py - zero,2)));
        dz = std::sqrt(arma::sum(arma::pow(pz - zero,2)));
    }




public:

    arma::cube     P,Ptmp,L;            // P(x, y, z) = R (probability)
    arma::mat      points; // (x, y, z) coordinates of each point
    arma::colvec3  x_ref_;          // point of reference of the grid
    BBox           bbox;



private:

    const pf::likelihood_model& likelihood_function;
    const pf::Measurement_h&    measurement_h;

    // position and orientation of the frame
    arma::colvec3  x_offset_;
    arma::colvec3  x_closet_ref_;

    std::size_t    N1,No;
    std::size_t    total_num_points; /// total number of grid points
    double         num_non_zero_points;  /// number of non-zero grid points

    GaussianKernel gauss_kernel;

    delta          delta_;
    length         length_;
    double         max_w;
    double         mean_w, sum_w;
    arma::mat33    Rot;
  //  arma::mat      hY;
    std::size_t    Y_dim;
    std::size_t    m_,n_,k_;
    std::chrono::time_point<std::chrono::steady_clock>          start_time, start_time2;
    std::chrono::time_point<std::chrono::steady_clock>          now_time;
    std::chrono::duration<double>                               diff_time;
    double                                                      update_dt, update_dt2;
    double                                                      percentage_non_zero;
    arma::colvec3                                               velocity;
    arma::colvec3                                               distance_travelled;
    double vel_conv_alpha;
    double bb_volume,bb_current_volume;
    double val_prior_conv, val_post_conv;
    double eps_trunc, eps_threashold;
    double delta_P;
    bool   bFirstSense;

    arma::rowvec3 w_mean;
    arma::mat33   w_cov;


    bool bTest;
    bool bUpdateColor;
    arma::colvec3   center;
    bool bYupdate;
    std::size_t counts;




};

}



#endif
