cmake_minimum_required(VERSION 3.22)
project(OpenTrustRegion VERSION 2.0.0 LANGUAGES Fortran C)

# get commands
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(otr OpenTrustRegion)  # Namespace

# ====  Options  ================================================================

# there are paths that work for both windows and unix but they repeat ${otr} :-(
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    set(${otr}_INSTALL_CMAKEDIR "${otr}/CMake"
        CACHE STRING "Directory to which CMake files are installed")
else()
    set(${otr}_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${otr}"
        CACHE STRING "Directory to which CMake files are installed")
endif()
message(STATUS "Showing option ${otr}_INSTALL_CMAKEDIR: ${${otr}_INSTALL_CMAKEDIR}")

option(${otr}_BUILD_TESTING "Build test-suite (was BUILD_TESTS)" ${PROJECT_IS_TOP_LEVEL})
message(STATUS "Showing option ${otr}_BUILD_TESTING: ${${otr}_BUILD_TESTING}")

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
message(STATUS "Showing option BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")

# ====  Build  ==================================================================


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
    if(BLAS_LIBRARIES)
        message(STATUS "Using user-specified BLAS")
    else()
        find_package(BLAS REQUIRED)
    endif()
    if(LAPACK_LIBRARIES)
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
if(BUILD_SHARED_LIBS)
    message(STATUS "Building OpenTrustRegion as a SHARED library")
else()
    message(STATUS "Building OpenTrustRegion as a STATIC library (default)")
endif()

# define source files
set(OPENTRUSTREGION_SOURCES src/opentrustregion.f90 src/c_interface.f90)

# add library
add_library(opentrustregion ${OPENTRUSTREGION_SOURCES})
add_library(${otr}::opentrustregion ALIAS opentrustregion)

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

if(${BUILD_SHARED_LIBS} AND APPLE)
    set_target_properties(opentrustregion PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

set(export_properties
  "${otr}_VERSION"
  )
set_target_properties(opentrustregion PROPERTIES ${otr}_VERSION ${${otr}_VERSION})
set_property(TARGET opentrustregion APPEND PROPERTY EXPORT_PROPERTIES "${export_properties}")

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
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# option to build testsuite
if (${otr}_BUILD_TESTING)
    # define test files
    set(OPENTRUSTREGION_TESTS
        tests/test_reference.f90
        tests/opentrustregion_unit_tests.f90
        tests/c_interface_unit_tests.f90
        tests/opentrustregion_mock.f90
        tests/c_interface_mock.f90
        tests/opentrustregion_system_tests.f90)

        # create shared test library which is then called from Python for testing
        add_library(testsuite SHARED ${OPENTRUSTREGION_TESTS})

        # link against BLAS and LAPACK and against OpenTrustRegion
        target_link_libraries(testsuite
                              PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
        )
endif()

# create targets file in the build tree for build-tree usage
export(TARGETS opentrustregion
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${otr}Targets.cmake"
       NAMESPACE "${otr}::"
)

# generate OpenTrustRegionConfig.cmake
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${otr}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${otr}Config.cmake
    INSTALL_DESTINATION ${${otr}_INSTALL_CMAKEDIR}
)

# generate version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${otr}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# install library
install(TARGETS opentrustregion
        EXPORT library_interface
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # Windows (.dll)
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # static libraries (.a)
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# install targets
install(EXPORT library_interface
        FILE ${otr}Targets.cmake
        NAMESPACE "${otr}::"
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
)

# install config files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}ConfigVersion.cmake
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
)
