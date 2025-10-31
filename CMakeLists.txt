cmake_minimum_required(VERSION 3.22)
project(OpenTrustRegion VERSION 2.0.0 LANGUAGES Fortran C)

# get commands
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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

# set integer size and define library suffix
message(STATUS "Using integer size ${BLA_SIZEOF_INTEGER}")
if(BLA_SIZEOF_INTEGER EQUAL 8)
    add_definitions(-DUSE_ILP64)
    set(LIB_SUFFIX "_64")
else()
    set(LIB_SUFFIX "_32")
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

# set properties
set_target_properties(
    opentrustregion PROPERTIES 
    # output name
    OUTPUT_NAME "opentrustregion${LIB_SUFFIX}" 
    # library version
    VERSION ${PROJECT_VERSION} 
    # interface version
    SOVERSION ${PROJECT_VERSION_MAJOR}
    # position-independent code - not automatically enabled for static libraries, but 
    # required when linking a static library into a shared library
    POSITION_INDEPENDENT_CODE ON 
)

# enable preprocessing
set_source_files_properties(${OPENTRUSTREGION_SOURCES} 
                            PROPERTIES Fortran_PREPROCESS ON
)

# link against BLAS and LAPACK
target_link_libraries(opentrustregion 
                      PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

# set include directories
target_include_directories(opentrustregion PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<INSTALL_INTERFACE:include>
)

# option to build testsuite
option(BUILD_TESTS "Build the testsuite" ON)

if (BUILD_TESTS)
    # create shared test library which is then called from Python for testing
    add_library(testsuite SHARED
                tests/test_reference.f90
                tests/opentrustregion_unit_tests.f90
                tests/c_interface_unit_tests.f90
                tests/opentrustregion_mock.f90
                tests/c_interface_mock.f90
                tests/opentrustregion_system_tests.f90
    )

    # link against BLAS and LAPACK and against OpenTrustRegion
    target_link_libraries(testsuite
                          PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
    )
endif()

# create targets file in the build tree for build-tree usage
export(TARGETS opentrustregion
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/OpenTrustRegionTargets.cmake"
       NAMESPACE OpenTrustRegion::
)

# generate OpenTrustRegionConfig.cmake
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/OpenTrustRegionConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenTrustRegion
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# generate version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/OpenTrustRegionConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# install library
install(TARGETS opentrustregion
        EXPORT OpenTrustRegionTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # static libraries (.a)
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ 
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install targets
install(EXPORT OpenTrustRegionTargets
        FILE OpenTrustRegionTargets.cmake
        NAMESPACE OpenTrustRegion::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenTrustRegion
)

# install config files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/OpenTrustRegionConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/OpenTrustRegionConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/OpenTrustRegion
)
