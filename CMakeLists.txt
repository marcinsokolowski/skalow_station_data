cmake_minimum_required(VERSION 3.5.1)
project(skalow_station_data VERSION 1.0.1)
include(GNUInstallDirs)

# file(GLOB skalow_station_data_sources "src/*.cpp")
# file(GLOB skalow_station_data_headers "src/*.h")
# add_library(skalow_station_data SHARED ${skalow_station_data_sources})

add_definitions(-D_UNIX -D_ENABLE_PROFILER_ -D_ENABLE_TRACING_ -D_ENABLE_VISUALIZATION_ -D_MONTE_CARLO_ -D_NO_ROOT_ -I${CMAKE_SOURCE_DIR}/src)

#set_target_properties(skalow_station_data PROPERTIES
#    VERSION ${PROJECT_VERSION}
#    SOVERSION 1
#    PUBLIC_HEADER "${skalow_station_data_headers}")
 
configure_file(skalow_station_data.pc.in skalow_station_data.pc @ONLY)


#install(TARGETS skalow_station_data
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#)
install(FILES ${CMAKE_BINARY_DIR}/skalow_station_data.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


find_library(CFITSIO_LIB cfitsio HINTS ENV LD_LIBRARY_PATH)
find_library(FFTW3_LIB fftw3 HINTS ENV FFTW_LIB REQUIRED)
find_library(LIBNOVA_LIB nova PATHS ENV LD_LIBRARY_PATH)
find_library(MSFITSLIB_LIB msfitslib PATHS ENV LD_LIBRARY_PATH)

# required to properly link the dynamic library :
# target_link_libraries(skalow_station_data ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(skalow_spectrometer apps/skalow_spectrometer.cpp)
target_link_libraries(skalow_spectrometer ${MSFITSLIB_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)


install(TARGETS skalow_spectrometer RUNTIME DESTINATION bin)

