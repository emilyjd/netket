


add_library(Sampler SHARED sampler.cc)
target_include_directories (Sampler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Sampler Machines Graph Hilbert Hamiltonian ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
