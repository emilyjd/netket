


add_library(Stepper SHARED stepper.cc)
target_include_directories (Stepper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Stepper Machines Observable Graph Hilbert Hamiltonian ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )

add_library(Learning SHARED learning.cc)
target_include_directories (Learning PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Learning Sampler Stepper Machines Observable Graph Hilbert Hamiltonian ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
