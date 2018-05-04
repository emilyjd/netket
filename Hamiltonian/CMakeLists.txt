


add_library(Hamiltonian SHARED hamiltonian.cc)

target_include_directories (Hamiltonian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Hamiltonian Hilbert Graph ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
