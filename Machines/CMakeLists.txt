


add_library(Machines SHARED machine.cc)
target_include_directories (Machines PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Machines Graph Hilbert Hamiltonian ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
