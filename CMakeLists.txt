
cmake_minimum_required(VERSION 2.8.3)
project(particle_filter)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  optitrack_rviz
  statistics_ml
  std_msgs
  message_generation
)

find_package(Armadillo)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -O3 -std=c++11 -mtune=native -fPIC -Wall -Wno-unused-variable")        ## Optimize
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

add_service_files(
  FILES
  String_cmd.srv
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES particle_filter
  CATKIN_DEPENDS roscpp
  DEPENDS system_lib
)


###########
## Build ##
###########

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

include_directories(${INCLUDE_DIR} ${catkin_INCLUDE_DIRS} ${INCLUDE_STATISTICS})


set(INC_FILES ${INCLUDE_DIR}/particle_filter/base_particle_filter.h
              ${INCLUDE_DIR}/particle_filter/particle_filter.h
              ${INCLUDE_DIR}/particle_filter/particle_filter_gmm.h
              ${INCLUDE_DIR}/particle_filter/sampling.h
              ${INCLUDE_DIR}/particle_filter/pf_visualisation.h
              ${INCLUDE_DIR}/particle_filter/particle_filter_definitions.h
              ${INCLUDE_DIR}/particle_filter/static_grid_filter.h
              ${INCLUDE_DIR}/particle_filter/reguliser.h
              ${INCLUDE_DIR}/point_mass_filter/point_mass_filter.h

)

set(SRC_FILES ${SRC_DIR}/base_particle_filter.cpp
              ${SRC_DIR}/particle_filter.cpp
              ${SRC_DIR}/particle_filter_gmm.cpp
              ${SRC_DIR}/sampling.cpp
              ${SRC_DIR}/static_grid_filter.cpp
              ${SRC_DIR}/reguliser.cpp
              ${SRC_DIR}/point_mass_filter/point_mass_filter.cpp
)


## Declare a cpp library
add_library(particle_filter
  ${SRC_FILES} ${INC_FILES}
)

target_link_libraries(particle_filter ${catkin_LIBRARIES} armadillo mlpack )

## Declare a cpp executable
add_executable(particle_filter_node
    ${SRC_DIR}/particle_filter_node.cpp
)

target_link_libraries(particle_filter_node
    particle_filter
)

## Add cmake target dependencies of the executable/library
## as an example, message headers may need to be generated before nodes
#add_dependencies(particle_filter particle_filter_generate_messages_cpp)

## Specify libraries to link a library or executable target against


#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS particle_filter particle_filter_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_particle_filter.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
