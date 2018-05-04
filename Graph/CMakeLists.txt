


add_library(Graph SHARED graph.cc)
target_include_directories (Graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Graph Hilbert ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
