#include "particle_filter/base_particle_filter.h"

namespace pf{

void  Base_particle_filter::set_visualisation_mode(opti_rviz::Vis_point_cloud::display_mode mode){
    vis_pf->set_display_type(mode);
}

void Base_particle_filter::set_color_mode(pf::color_type color_t){
    this->color_t = color_t;
}

}
