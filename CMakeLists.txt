cmake_minimum_required(VERSION 3.22)
project(OpenTrustRegion LANGUAGES Fortran C)
include(GNUInstallDirs)

# default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# set flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
endif()

# user-supplied integer size
set(INTEGER_SIZE "" CACHE STRING "Size of integers for library (4=32-bit, 8=64-bit)")

# user-supplied BLAS and LAPACK libraries
set(BLAS_LIBRARIES "" CACHE STRING "User-specified BLAS library path(s)")
set(LAPACK_LIBRARIES "" CACHE STRING "User-specified LAPACK library path(s)")

# check if integer size is provided
if(INTEGER_SIZE)
    set(BLA_SIZEOF_INTEGER ${INTEGER_SIZE} CACHE INTERNAL "BLAS/LAPACK integer size")
    if(BLAS_LIBARIES)
        message(STATUS "Using user-specified BLAS")
    else()
        find_package(BLAS REQUIRED)
    endif()
    if(LAPACK_LIBARIES)
        message(STATUS "Using user-specified LAPACK")
    else()
        find_package(LAPACK REQUIRED)
    endif()
else()
    if(BLAS_LIBRARIES OR LAPACK_LIBRARIES)
        message(FATAL_ERROR
            "You provided explicit BLAS or LAPACK libraries but did not specify INTEGER_SIZE. "
            "Please set INTEGER_SIZE to 4 (LP64) or 8 (ILP64) to ensure compatibility."
        )
    else()
        # autodetect BLAS and LAPACK, try 32-bit first, then 64-bit
        set(BLA_SIZEOF_INTEGER 4 CACHE INTERNAL "BLAS/LAPACK integer size")
        find_package(BLAS)
        find_package(LAPACK)
        if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
            set(BLA_SIZEOF_INTEGER 8 CACHE INTERNAL "BLAS/LAPACK integer size")
            find_package(BLAS REQUIRED)
            find_package(LAPACK REQUIRED)
        endif()
    endif()
endif()

message(STATUS "Using integer size ${BLA_SIZEOF_INTEGER}")
if(BLA_SIZEOF_INTEGER EQUAL 8)
    add_definitions(-DUSE_ILP64)
    set(LIB_SUFFIX "64")
else()
    set(LIB_SUFFIX "32")
endif()

# default to static libraries
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
if(BUILD_SHARED_LIBS)
    message(STATUS "Building OpenTrustRegion as a SHARED library")
else()
    message(STATUS "Building OpenTrustRegion as a STATIC library (default)")
endif()

# define source files
set(OPENTRUSTREGION_SOURCES src/opentrustregion.f90 src/c_interface.f90)

# add library
add_library(opentrustregion ${OPENTRUSTREGION_SOURCES})

# set library version
set_target_properties(opentrustregion PROPERTIES VERSION 2.0.0 SOVERSION 2)

# compile position-independent code - this is not automatically enabled for static 
# libraries, but is required when linking a static library into a shared library
set_target_properties(opentrustregion PROPERTIES
    OUTPUT_NAME "opentrustregion${LIB_SUFFIX}"
    POSITION_INDEPENDENT_CODE ON
)

# enable preprocessing
set_source_files_properties(${OPENTRUSTREGION_SOURCES} PROPERTIES Fortran_PREPROCESS ON)

# link against BLAS and LAPACK
target_link_libraries(opentrustregion PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# set include directories
target_include_directories(opentrustregion PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# option to build testsuite
option(BUILD_TESTS "Build the testsuite" ON)

if (BUILD_TESTS)
    # create shared test library which is then called from Python for testing
    add_library(
        testsuite SHARED
        tests/opentrustregion_unit_tests.f90
        tests/c_interface_unit_tests.f90
        tests/opentrustregion_mock.f90
        tests/c_interface_mock.f90
        tests/opentrustregion_system_tests.f90
    )
  
    # link against BLAS and LAPACK and against OpenTrustRegion
    target_link_libraries(
        testsuite
        PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
    )
    
    # enable preprocessing
    set_source_files_properties(tests/opentrustregion_system_tests.f90 PROPERTIES Fortran_PREPROCESS ON)
endif()

# set paths to find module files, header files and library in build directory
set(OpenTrustRegion_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/ ${CMAKE_BINARY_DIR})

set(OpenTrustRegion_LIBRARY $<TARGET_FILE:opentrustregion>)

# create the OpenTrustRegionConfig.cmake from the template
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/OpenTrustRegionConfig.cmake
  @ONLY
)

# export the target to make OpenTrustRegion::OpenTrustRegion usable
install(TARGETS opentrustregion
        EXPORT OpenTrustRegionTargets
        LIBRARY DESTINATION lib # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION lib # static libraries (.a)
        INCLUDES DESTINATION include)
install(FILES include/opentrustregion.h DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/)

# export the target with namespace OpenTrustRegion::
install(EXPORT OpenTrustRegionTargets
        NAMESPACE OpenTrustRegion::
        DESTINATION lib/cmake/OpenTrustRegion)

# generate and install the config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/OpenTrustRegionConfig.cmake
    INSTALL_DESTINATION lib/cmake/OpenTrustRegion
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OpenTrustRegionConfig.cmake
        DESTINATION lib/cmake/OpenTrustRegion)
