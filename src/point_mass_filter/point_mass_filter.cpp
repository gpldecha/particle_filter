#include "point_mass_filter/point_mass_filter.h"
#include <tf/transform_broadcaster.h>
#include <algorithms/convolution.h>
#include <assert.h>

namespace pf{

int throt_time = 2;


Point_mass_filter::Point_mass_filter(const pf::likelihood_model& likelihood_function,
                                     const Measurement_h &measurement_h,
                                     delta delta_,length length_,std::size_t Y_dim)
    :likelihood_function(likelihood_function),
      measurement_h(measurement_h),
      Y_dim(Y_dim)
{

    arma::colvec3 center;
    center.zeros();
    reset(center,delta_,length_);

}

void Point_mass_filter::reset(const arma::colvec3& center,delta& delta_,length& length_){

    this->center  = center;
    this->delta_  = delta_;
    this->length_ = length_;

    bTest = true;
    bFirstSense  = false;
    bUpdateColor = true;
    bYupdate     = true;
    counts       = 0;

    // number of grid points x-axis
    m_ = std::floor(length_.m/delta_.m);
    // number of grid points y-axis
    n_ = std::floor(length_.n/delta_.n);
    // number of grid points z-axis
    k_ = std::floor(length_.k/delta_.k);

    if(m_ % 2 == 0){m_ = m_ + 1;}
    if(n_ % 2 == 0){n_ = n_ + 1;}
    if(k_ % 2 == 0){k_ = k_ + 1;}

    x_offset_(0) = (0.5 * static_cast<double>(m_-1) * delta_.m);
    x_offset_(1) = (0.5 * static_cast<double>(n_-1) * delta_.n);
    x_offset_(2) = (0.5 * static_cast<double>(k_-1) * delta_.k);

    x_ref_.zeros();
    x_ref_ = center -  x_offset_;

    P = arma::cube(m_,n_,k_);
    P = P.ones();
    P = P / static_cast<double>(m_ * n_ * k_);

    L.resize(m_,n_,k_);

    gauss_kernel.var_x = std::pow(0.05,2.0);
    gauss_kernel.var_y = std::pow(0.05,2.0);
    gauss_kernel.var_z = std::pow(0.05,2.0);

    gauss_kernel.reset(delta_);
    gauss_kernel.print();


    colors.resize(P.n_elem);
   // hY.resize(P.n_elem,Y_dim);


    points = arma::zeros<arma::mat>(P.n_elem,3);
    max_w  = P.max();
    get_coordiantes();

    arma::uword index;
    arma::sqrt(arma::sum(   arma::pow(points - arma::repmat(center.st(),points.n_rows,1),2) ,1)   ).min(index);
    x_closet_ref_ = points.row(index).st();

    x_ref_ = x_ref_ + (center - x_closet_ref_);

    get_coordiantes();
    arma::sqrt(arma::sum(   arma::pow(points - arma::repmat(center.st(),points.n_rows,1),2) ,1)   ).min(index);
    x_closet_ref_ = points.row(index).st();

    update_dt  = 0.1;
    update_dt2 = 0.05;
    eps_trunc  = 0.05;

    No = 500;
    N1 = 100000;

    get_bbox();
    bb_current_volume = bb_volume;


}

void Point_mass_filter::update(const arma::colvec &u, const arma::colvec &Y, double duration){
    bool debug     = false;
    bool debug_lik = false;

    if(debug_lik)
    {
        update_lik_debug(u,Y,duration);
    }else{
        if(debug)
        {
            update_debug(u,Y,duration);
        }else{
            update_real(u,Y,duration);
        }
    }
}

void Point_mass_filter::update_debug(const arma::colvec &u, const arma::colvec &Y, double duration){
    now_time   = std::chrono::steady_clock::now();
    arma::colvec Y_ = Y;
    Y_(1) = 0;
    Y_(2) = 0;

    diff_time = now_time-start_time;
    if(diff_time.count() > 0.001)
    {

        gauss_kernel.var_x = std::pow(0.001,2);
        gauss_kernel.var_y = std::pow(0.1,2);
        gauss_kernel.var_z = std::pow(0.1,2);
        gauss_kernel.reset(delta_);
        gauss_kernel.print();

        assert(gauss_kernel.kernel_x.n_elem >= 3);
        assert(gauss_kernel.kernel_y.n_elem >= 3);
        assert(gauss_kernel.kernel_z.n_elem >= 3);


        total_num_points = P.n_elem;
        eps_threashold = eps_trunc / static_cast<double>(total_num_points);
        assert(std::isfinite(eps_threashold));
        get_num_non_zero();

        std::cout<< "eps_threashold: " << eps_threashold << std::endl;


        normalise();
        /*   Ptmp = P;
         stats::conv3<double>(P,Ptmp,gauss_kernel.kernel_x,gauss_kernel.kernel_y,gauss_kernel.kernel_z);

         double max_r,max_c, max_k;
         get_boundary_values(max_r,max_c,max_k);

        //if(max_r > 0.1){
        //    increase_lenght(2,0,0);
       // }
        if(max_c > 0.1){
            increase_lenght(0,2,0);
        }
        if(max_k > 0.1){
            increase_lenght(0,0,2);
        }
*/
        normalise();
        get_num_non_zero();
        std::cout<< "num_non_zero_points: " << num_non_zero_points << std::endl;
        get_bbox();

        std::cout<< "bb_volume:         " << bb_volume << std::endl;
        std::cout<< "bb_current_volume: " << bb_current_volume << std::endl;
        std::cout<< "volume:            " << bb_volume/bb_current_volume << std::endl;

        assert(std::isfinite(bb_volume/bb_current_volume));

        if(bb_volume/bb_current_volume < 0.8)
        {
            std::cout<< "start TRANSFORM BBOX" << std::endl;
            transform_bbox();
          //  hY.resize(P.n_elem,Y_dim);
            get_coordiantes();
            //get_num_non_zero();
            std::cout<< "end TRANSFORM BBOX" << std::endl;
        }
        /*
       int num_if_increase = get_num_if_increase();
       if((num_if_increase < N1) && (num_non_zero_points < No))
       {
           ROS_INFO_STREAM("increase_densit()");
           increase_density();
           hY.resize(P.n_elem,3);
       }else
       {
            /// Check if need to decrease density of points
           if( (num_non_zero_points > N1)){
               ROS_INFO_STREAM("decrease_density()");
               decrease_density();
               hY.resize(P.n_elem,3);
           }
       }

       ROS_INFO_STREAM_THROTTLE(throt_time,"Nt: " << total_num_points);*/
        start_time = std::chrono::steady_clock::now();
        max_w = P.max();
    }

}

void Point_mass_filter::update_lik_debug(const arma::colvec &u, const arma::colvec &Y, double duration){
    arma::colvec Y_ = Y;

     if(Y_(0) != 0){
        bFirstSense=true;
    }

    //motion_update(u);

    if(bYupdate){
        measurement_update(Y_);
    }
}

void Point_mass_filter::update_real(const arma::colvec& u, const arma::colvec& Y, double duration){
    now_time   = std::chrono::steady_clock::now();

    arma::colvec Y_ = Y;

    if(Y_(0) != 0){
        bFirstSense=true;
    }

    motion_update(u);

    if(bYupdate){

        diff_time = now_time - start_time2;
        if(diff_time.count() > update_dt2){
            measurement_update(Y_);
            start_time2 = std::chrono::steady_clock::now();
        }
    }

    diff_time = now_time-start_time;
    if(diff_time.count() > update_dt)
    {

        get_num_non_zero();
        get_bbox();

       // ROS_INFO_STREAM_THROTTLE(throt_time,"    ");
      //  ROS_INFO_STREAM_THROTTLE(throt_time,"Nt: " << total_num_points);
      //  ROS_INFO_STREAM_THROTTLE(throt_time,"N:  " << num_non_zero_points);
     // ROS_INFO_STREAM_THROTTLE(throt_time,"bb_volume:  " << bb_volume);
      //  ROS_INFO_STREAM_THROTTLE(throt_time,"%vol:  " << bb_volume/bb_current_volume);

    /*    if(bb_volume/bb_current_volume < 0.4)
        {
            std::cout<< "start TRANSFORM BBOX" << std::endl;
            transform_bbox();
            hY.resize(P.n_elem,Y_dim);
            get_coordiantes();
            std::cout<< "end TRANSFORM BBOX" << std::endl;
        }

        int num_if_increase = get_num_if_increase();
        if((num_if_increase < N1) && (num_non_zero_points < No))
        {
          //  ROS_INFO_STREAM("increase_densit()");
            increase_density();
            hY.resize(P.n_elem,Y_dim);
        }else
        {
            /// Check if need to decrease density of points
            if( (num_non_zero_points > N1)){
              //  ROS_INFO_STREAM("decrease_density()");
                decrease_density();
                hY.resize(P.n_elem,Y_dim);
            }
        }*/

        start_time = std::chrono::steady_clock::now();
        bYupdate = true;
    }

   // ROS_INFO_STREAM_THROTTLE(throt_time,"Nt: " << total_num_points);
   // ROS_INFO_STREAM_THROTTLE(throt_time,"N : " << num_non_zero_points);
}

void Point_mass_filter::motion_update(const arma::colvec &u, double duration){

    x_ref_(0) = x_ref_(0) + u(0);
    x_ref_(1) = x_ref_(1) + u(1);
    x_ref_(2) = x_ref_(2) + u(2);

    distance_travelled = distance_travelled + u;

    if(bFirstSense && arma::norm(distance_travelled) >= 0.03)
    {
        diff_time   = now_time-start_time2;

        velocity = distance_travelled / diff_time.count();
        vel_conv_alpha = 0.001;
        gauss_kernel.var_x = vel_conv_alpha * std::pow(velocity(0),2);
        gauss_kernel.var_y = vel_conv_alpha * std::pow(velocity(1),2);
        gauss_kernel.var_z = vel_conv_alpha * std::pow(velocity(2),2);
        gauss_kernel.reset(delta_);
        gauss_kernel.print();
        assert(gauss_kernel.kernel_x.n_elem >= 3);
        assert(gauss_kernel.kernel_y.n_elem >= 3);
        assert(gauss_kernel.kernel_z.n_elem >= 3);

        get_num_non_zero();
        Ptmp = P;
        stats::conv3<double>(P,Ptmp,gauss_kernel.kernel_x,gauss_kernel.kernel_y,gauss_kernel.kernel_z);
        normalise();
        double max_r,max_c, max_k;
        get_boundary_values(max_r,max_c,max_k);
        if(max_c > 0.1 || P.n_cols < 3){

            increase_lenght(0,2,0);
        }
        if(max_k > 0.1 || P.n_cols < 3){
            increase_lenght(0,0,2);
        }

        normalise();
        distance_travelled.zeros();
        start_time2 = std::chrono::steady_clock::now();
    }
}

void Point_mass_filter::measurement_update(const arma::colvec& Y){



    get_coordiantes();
    //measurement_h(hY,points,Rot); // bottelneck

    if(L.n_elem != P.n_elem)
    {
        L.ones(P.n_rows,P.n_cols,P.n_slices);
    }else{
        L.ones();
    }

    likelihood_function(L.memptr(),Y,points,Rot);
    P = L % P;
    normalise();

}

void Point_mass_filter::normalise(){
    double sum_denom = 0;
    double delta_cube = (delta_.m * delta_.n * delta_.k);
    P = P * delta_cube;
    sum_denom = arma::sum(arma::vectorise(P));
    if(sum_denom != 0)
    {
        P = P / sum_denom;
    }else{
        std::cout<< "pmf normalise() error sum denom == 0" << std::endl;
        P.ones();
        normalise();
    }
    max_w = arma::max(arma::vectorise(P));
}

void Point_mass_filter::interpolate(){
    Ptmp = P;

    int rows   = static_cast<int>(P.n_rows);
    int cols   = static_cast<int>(P.n_cols);
    int slices = static_cast<int>(P.n_slices);

    // X-axis interpolation
    if(rows >= 3){
        for(int c = 0; c < cols; c++ ){
            for(int k = 0; k < slices; k++){

                P(0,c,k)      = 0.5 * Ptmp(0,c,k)       + 0.5 * Ptmp(1,c,k);
                P(rows-1,c,k) = 0.5 * Ptmp(rows-2,c,k)  + 0.5 * Ptmp(rows-1,c,k);
                for(int r = 1; r < rows-1;r++){
                    P(r,c,k)  = 0.5 * Ptmp(r-1,c,k) + Ptmp(r,c,k)  + 0.5 * Ptmp(r+1,c,k);
                }
            }
        }
        Ptmp = P;
    }
    // Y-axis interpolation
    if(cols >= 3){
        for(int r = 0; r < rows;r++){
            for(int k = 0; k < slices; k++){

                P(r,0,k)      = 0.5 * Ptmp(r,0,k)       + 0.5 * Ptmp(r,1,k);
                P(r,cols-1,k) = 0.5 * Ptmp(r,cols-2,k)  + 0.5 * Ptmp(r,cols-1,k);
                for(int c = 1; c < cols-1; c++ ){
                    P(r,c,k) = 0.5 * Ptmp(r,c-1,k) + Ptmp(r,c,k)  + 0.5 * Ptmp(r,c+1,k);
                }
            }
        }
        Ptmp = P;
    }
    // Z-axis interpolation
    if(slices >= 3){
        for(int r = 0; r < rows;r++){
            for(int c = 0; c < cols; c++ ){
                P(r,c,0)        = 0.5 * Ptmp(r,c,0)         + 0.5 * Ptmp(r,c,1);
                P(r,c,slices-1) = 0.5 * Ptmp(r,c,slices-2)  + 0.5 * Ptmp(r,c,slices-1);
                for(int k = 1; k < slices-1; k++){
                    P(r,c,k) = 0.5 * Ptmp(r,c,k-1) + Ptmp(r,c,k)  + 0.5 * Ptmp(r,c,k+1);
                }
            }
        }
    }
}

void Point_mass_filter::increase_density(){

    // interpolata
    Ptmp.zeros(2*m_-1,2*n_-1,2*k_-1);

    arma::colvec3 point_zero;
    arma::colvec3 point_zero2;
    point_zero = points.row(0).st();

    std::size_t i, j, k;
    int n_rows = static_cast<int>(P.n_rows);
    int n_cols = static_cast<int>(P.n_cols);
    int RC = n_rows * n_cols;
    int sub;

    for(int n = 0; n < static_cast<int>(P.n_elem);n++){
        // convert linear index to r, c, s
        k   = n / RC; // slice value (Note: integer division)
        sub = n % RC;
        j = sub / n_rows; // column value (Note: integer division)
        i = sub % n_rows; // row value

        Ptmp(2*i,2*j,2*k) = P.at(n);
    }

    P         = Ptmp;
    points.resize(P.n_elem,3);
    //hY.resize(P.n_elem,Y_dim);
    m_       = P.n_rows;
    n_       = P.n_cols;
    k_       = P.n_slices;

    delta_.m = delta_.m/2;
    delta_.n = delta_.n/2;
    delta_.k = delta_.k/2;
    colors.resize(P.n_elem);
    vis_pf->initialise("world",points);

    get_coordiantes();

    point_zero2(0) = delta_.m + x_ref_(0);
    point_zero2(1) = delta_.n + x_ref_(1);
    point_zero2(2) = delta_.k + x_ref_(2);


    x_ref_ = x_ref_ + (point_zero - point_zero2);

    interpolate();

    normalise();
}

void Point_mass_filter::decrease_density(){


    m_ = P.n_rows;
    n_ = P.n_cols;
    k_ = P.n_slices;

    std::size_t m_n = (m_ - 1.0)/2.0 + 1;
    std::size_t n_n = (n_ - 1.0)/2.0 + 1;
    std::size_t k_n = (k_ - 1.0)/2.0 + 1;
    if(m_n <= 2){m_n = 1;}
    if(n_n <= 2){n_n = 1;}
    if(k_n <= 2){k_n = 1;}

    std::cout<< "m:   " << m_ << " " << n_ << " " << k_ << std::endl;
    std::cout<< "m_n: " << m_n << " " << n_n << " " <<  k_n<< std::endl;


    Ptmp.zeros(m_n,n_n,k_n);

    arma::colvec3 point_zero;
    arma::colvec3 point_zero2;

    // P(1,1,1) -> P(0,0,0)
    point_zero(0) = delta_.m + x_ref_(0);
    point_zero(1) = delta_.n + x_ref_(1);
    point_zero(2) = delta_.k + x_ref_(2);
    // point_zero = points.row(0).st();

    std::size_t i, j, k;
    int n_rows = static_cast<int>(Ptmp.n_rows);
    int n_cols = static_cast<int>(Ptmp.n_cols);
    int RC = n_rows * n_cols;
    int sub;
    std::cout<< "before for loop" << std::endl;
    for(int n = 0; n < static_cast<int>(Ptmp.n_elem);n++){
        // convert linear index to r, c, s
        k   = n / RC; // slice value (Note: integer division)
        sub = n % RC;
        j = sub / n_rows; // column value (Note: integer division)
        i = sub % n_rows; // row value
        i = 2*i;
        j = 2*j;
        k = 2*k;
        Ptmp.at(n) = P(i,j,k);
    }
    P = Ptmp;
    points.set_size(P.n_elem,3);
    m_ = P.n_rows;
    n_ = P.n_cols;
    k_ = P.n_slices;

    delta_.m = 2 * delta_.m;
    delta_.n = 2 * delta_.n;
    delta_.k = 2 * delta_.k;

    colors.resize(P.n_elem);
    vis_pf->initialise("world",points);

    point_zero2(0) = delta_.m + x_ref_(0);
    point_zero2(1) = delta_.n + x_ref_(1);
    point_zero2(2) = delta_.k + x_ref_(2);
    // apply correction
    x_ref_ = x_ref_ + (point_zero2 - point_zero);

    // std::cout<< "P : (" << P.n_rows << " x " << P.n_cols << " x " << P.n_slices << ") " << std::endl;


    std::cout<< "after loop" << std::endl;


}

void Point_mass_filter::increase_lenght(int dm, int dn, int dk){

    Ptmp.zeros(m_ + dm,n_ + dn,k_ + dk);

    double delta_vol = delta_.k * delta_.m * delta_.n;
   // ROS_INFO_STREAM_THROTTLE(throt_time,"max_w: " << max_w);
    double eps       = 0.01 * max_w;
   // std::cout<< "eps: " << eps << std::endl;
     Ptmp = Ptmp * eps;


    std::size_t i, j, k;
    int n_rows = static_cast<int>(P.n_rows);
    int n_cols = static_cast<int>(P.n_cols);
    int RC = n_rows * n_cols;
    int sub;

    std::size_t ii,jj,kk;
    ii = dm/2;
    jj = dn/2;
    kk = dk/2;

    for(int n = 0; n < static_cast<int>(P.n_elem);n++){
        // convert linear index to r, c, s
        k   = n / RC; // slice value (Note: integer division)
        sub = n % RC;
        j = sub / n_rows; // column value (Note: integer division)
        i = sub % n_rows; // row value

        Ptmp(i+ii,j+jj,k+kk) = P.at(n);
    }

    P = Ptmp;

    // iterate over all boundary points



    points.resize(P.n_elem,3);
  //  hY.resize(P.n_elem,Y_dim);

    m_       = P.n_rows;
    n_       = P.n_cols;
    k_       = P.n_slices;
    colors.resize(P.n_elem);
    vis_pf->initialise("world",points);

    x_ref_(0) = x_ref_(0) - static_cast<double>(dm)  *  delta_.m/2.0;
    x_ref_(1) = x_ref_(1) - static_cast<double>(dn)  *  delta_.n/2.0;
    x_ref_(2) = x_ref_(2) - static_cast<double>(dk)  *  delta_.k/2.0;

    // has to be called after the x_ref_ was changed
    get_coordiantes();
}

void Point_mass_filter::get_bbox(){

    double max_x = -std::numeric_limits<double>::max();
    double min_x =  std::numeric_limits<double>::max();

    double max_y = max_x;
    double min_y = min_x;

    double max_z = max_x;
    double min_z = min_x;

    // assumed aligned with the axis
    for(std::size_t n = 0; n < points.n_rows;n++)
    {
        if(P.at(n) != 0){

            if(points(n,0) >= max_x)
            {
                bbox.max_x_n = n;
                max_x        = points(n,0);
            }
            if(points(n,0) <= min_x)
            {
                bbox.min_x_n = n;
                min_x        = points(n,0);
            }
            // find y max-min
            if(points(n,1) >= max_y)
            {
                bbox.max_y_n = n;
                max_y        = points(n,1);
            }
            if(points(n,1) <= min_y)
            {
                bbox.min_y_n = n;
                min_y        = points(n,1);
            }
            // find z max-min
            if(points(n,2) >= max_z)
            {
                bbox.max_z_n = n;
                max_z        = points(n,2);
            }
            if(points(n,2) <= min_z)
            {
                bbox.min_z_n = n;
                min_z        = points(n,2);
            }

        }
    }

    bbox.update(points,P,delta_);


    bb_volume = bbox.dist_x * bbox.dist_y * bbox.dist_z;

    arma::colvec3 c_bb = points.row(0).st();
    arma::colvec3 c_tt = points.row(points.n_rows-1).st();

    double c_dist_x = std::fabs(c_tt(0) - c_bb(0) + delta_.m);
    double c_dist_y = std::fabs(c_tt(1) - c_bb(1) + delta_.n);
    double c_dist_z = std::fabs(c_tt(2) - c_bb(2) + delta_.k);


    bb_current_volume = c_dist_x * c_dist_y * c_dist_z;

    assert(std::isfinite(bb_volume/bb_current_volume));

    tf_debuf(bbox.tt,"tt");
    tf_debuf(bbox.bb,"bb");

    tf_debuf(c_tt,"c_tt");
    tf_debuf(c_bb,"c_bb");

}

void Point_mass_filter::get_index_boundary(int& max, int& min,int n, int num_cells){

    if(n < 3){  // if the new bounding box gets less than 3 grid cells for this dimension
        if(n == 2)
        {
            if(max < num_cells- 2)
            {
                max = max + 1;
            }

        }else if(n == 1){

            // if space on bow sides add one to each
            if(max < num_cells-2 && min > 0)
            {
                max = max + 1;
                min = min - 1;
                // only space on min side
            }else if(max == num_cells-1 && min > 1)
            {
                min = min - 2;
                // only space on max side
            }else if(max < num_cells-3 && min == 0)
            {
                max = max + 2;
            }
        }
    }
}

void Point_mass_filter::transform_bbox(){

    std::size_t ii,jj,kk;
    ii = jj = kk = 0;

    int r_max = bbox.r_max;
    int r_min = bbox.r_min;
    int c_max = bbox.c_max;
    int c_min = bbox.c_min;
    int k_max = bbox.k_max;
    int k_min = bbox.k_min;

    int m_n = r_max - r_min + 1;
    int n_n = c_max - c_min + 1;
    int k_n = k_max - k_min + 1;

    get_index_boundary(r_max,r_min,m_n,P.n_rows);
    get_index_boundary(c_max,c_min,n_n,P.n_cols);
    get_index_boundary(k_max,k_min,k_n,P.n_slices);

    m_n = r_max - r_min + 1;
    n_n = c_max - c_min + 1;
    k_n = k_max - k_min + 1;

    std::cout<< "r: " << r_max << " " << r_min << std::endl;
    std::cout<< "r: " << c_max << " " << c_min << std::endl;
    std::cout<< "r: " << k_max << " " << k_min << std::endl;
    std::cout<< "n: " << m_n << " " << n_n << " " << k_n << std::endl;

    std::cout<< "Ptmp.reshape"<< std::endl;
    Ptmp.zeros(m_n,n_n,k_n);


    std::cout<< "m_ " << m_ << std::endl;
    std::cout<< "n_ " << n_ << std::endl;
    std::cout<< "k_ " << k_ << std::endl;
    std::cout<< "m_n: " << m_n << std::endl;
    std::cout<< "n_n: " << n_n << std::endl;
    std::cout<< "k_n: " << k_n << std::endl;

    for(int r = r_min; r <= r_max;r++){
        for(int c = c_min; c <= c_max; c++){
            for(int k =  k_min; k <= k_max;k++){

                //   std::cout<< "("<<ii<<","<<jj<<","<<kk<<") <---- " << "(" << r << "," << c << "," << k << ")" << std::endl;
                assert(ii >= 0 && ii < Ptmp.n_rows);
                assert(jj >= 0 && jj < Ptmp.n_cols);
                assert(kk >= 0 && kk < Ptmp.n_slices);
                Ptmp(ii,jj,kk) = P(r,c,k);
                kk = kk + 1;
            }
            jj = jj + 1;
            kk = 0;
        }
        jj = 0;
        ii = ii + 1;
    }

    P = Ptmp;

    m_ = P.n_rows;
    // number of grid points y-axis
    n_ = P.n_cols;
    // number of grid points z-axis
    k_ = P.n_slices;

    x_ref_ = bbox.bb;

    points.set_size(P.n_elem,3);
    get_coordiantes();

    colors.resize(P.n_elem);
    vis_pf->initialise("world",points);

}

void Point_mass_filter::get_coordiantes(){

    std::size_t i, j, k;
    int n_rows = static_cast<int>(P.n_rows);
    int n_cols = static_cast<int>(P.n_cols);
    int RC = n_rows * n_cols;
    int sub;
    for(int n = 0; n < static_cast<int>(P.n_elem);n++){
        // convert linear index to r, c, s
        k   = n / RC; // slice value (Note: integer division)
        sub = n % RC;

        j = sub / n_rows; // column value (Note: integer division)
        i = sub % n_rows; // row value

        points(n,0) = i * delta_.m + x_ref_(0);
        points(n,1) = j * delta_.n + x_ref_(1);
        points(n,2) = k * delta_.k + x_ref_(2);
    }
}

void Point_mass_filter::set_rotation(const arma::mat33& Rot){
    this->Rot = Rot;
}

void Point_mass_filter::print_prob() const {
    std::cout<< "=== P === " << std::endl;
    std::cout<< "( " << P.n_rows << " x " << P.n_cols << " x " << P.n_slices << ") " << std::endl;
    std::cout<< "m : " << m_ << " n: " << n_ << " k: " << k_ << std::endl;
}

double Point_mass_filter::mean_weight() const{
    return arma::sum(arma::vectorise(P)) / static_cast<double>(P.n_elem);
}

void Point_mass_filter::init_visualise(ros::NodeHandle& node,const std::string& topic_name){
    vis_pf.reset(new opti_rviz::Vis_point_cloud(node,topic_name));
  // vis_pf->set_channel(opti_rviz::Vis_point_cloud::CHANNEL_TYPE::Intensity);
    vis_pf->initialise("world",points);
    vis_pf->set_display_type(opti_rviz::Vis_point_cloud::ONLY_HIGH_WEIGHTS);
}

void Point_mass_filter::compute_color(pf::color_type c_type){
    if(c_type == C_LIKE){
        //  std::cout<< "c_like" << std::endl;
        double max_L = arma::max(arma::vectorise(L));
        for(std::size_t i = 0 ; i < L.n_elem;i++){
            ColorMap::jetColorMap(rgb,L.at(i),0,max_L);
            colors[i][0]    = ((float)rgb[0])/255;
            colors[i][1]    = ((float)rgb[1])/255;
            colors[i][2]    = ((float)rgb[2])/255;
        }
    }else{
        for(std::size_t i = 0 ; i < P.n_elem;i++){
            ColorMap::jetColorMap(rgb,P.at(i),0,max_w);
            colors[i][0]    = ((float)rgb[0])/255;
            colors[i][1]    = ((float)rgb[1])/255;
            colors[i][2]    = ((float)rgb[2])/255;
        }
    }
}

void Point_mass_filter::visualise(){
    //  std::cout<< "compute_color" << std::endl;
    compute_color(color_t);
    // std::cout<< "update" << std::endl;
    vis_pf->update(points,colors,P.memptr(),0.8 * max_w);
    vis_pf->publish();

    static tf::TransformBroadcaster br;
    tf::Transform transform;
    transform.setOrigin( tf::Vector3(x_closet_ref_(0),x_closet_ref_(1),x_closet_ref_(2)) );
    tf::Quaternion q;
    q.setRPY(0, 0, 0);
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "x_closet_ref_"));
}


void Point_mass_filter::get_boundary_values(double& max_r,double& max_c, double& max_k){

    max_r = 0;
    max_c = 0;
    max_k = 0;
    // get X-boundary values
    std::size_t r = 0;
    std::size_t c = 0;
    std::size_t d = 0;
    std::size_t r_end = P.n_rows-1;
    std::size_t c_end = P.n_cols-1;
    std::size_t d_end = P.n_slices-1;


    for(c=0; c < P.n_cols;c++){
        for(d=0; d < P.n_slices;d++){

            if(max_r < P(0,c,d))
            {
                max_r = P(0,c,d);
            }

            if(max_r < P(r_end,c,d))
            {
                max_r= P(r_end,c,d);
            }

        }
    }

    r = 0;
    for(r=0;r < P.n_rows;r++){
        // iterating over D
        for(c=0; c < P.n_cols;c++){

            if(max_k < P(r,c,0))
            {
                max_k = P(r,c,0);
            }
            if(max_k < P(r,c,d_end))
            {
                max_k = P(r,c,d_end);
            }
        }
        for(d = 0; d < P.n_slices;d++)
        {
            if(max_c < P(r,0,d))
            {
                max_c = P(r,0,d);
            }
            if(max_c < P(r,c_end,d))
            {
                max_c = P(r,c_end,d);
            }
        }
    }
    //  std::cout<< "max__r " << max_r << "  max: " << max_w << std::endl,

    max_r = max_r / max_w;
    max_c = max_c / max_w;
    max_k = max_k / max_w;

    std::cout<< "max_r: " <<  max_r << std::endl;
    std::cout<< "max_c: " <<  max_c << std::endl;
    std::cout<< "max_k: " <<  max_k << std::endl;


}


const Point_mass_filter::delta& Point_mass_filter::get_delta() const{
    return delta_;
}

}

