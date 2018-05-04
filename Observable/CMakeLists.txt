


add_library(Observable SHARED observable.cc)
target_include_directories (Observable PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Observable Hilbert ${MPI_LIBRARIES})
set( CMAKE_VERBOSE_MAKEFILE on )
