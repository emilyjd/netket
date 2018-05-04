


add_library(Hilbert SHARED hilbert.cc)
target_include_directories (Hilbert PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Hilbert ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
