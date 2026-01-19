cmake_minimum_required(VERSION 3.22)
project(OpenTrustRegion VERSION 2.0.0 LANGUAGES Fortran C)

# get commands
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckFortranFunctionExists)

# set namespace
set(otr OpenTrustRegion)

# set components
set(runtime_component runtime)
set(test_component testsuite)

# ====  Options  ================================================================

# set installation directories, this needs to be separate for windows and unix to avoid 
# repeating ${otr} twice in the path
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    set(${otr}_INSTALL_CMAKEDIR "${otr}/CMake"
        CACHE STRING "Directory to which CMake files are installed")
else()
    set(${otr}_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${otr}"
        CACHE STRING "Directory to which CMake files are installed")
endif()
message(STATUS "Showing option ${otr}_INSTALL_CMAKEDIR: ${${otr}_INSTALL_CMAKEDIR}")

# default to not building tests
option(${otr}_BUILD_TESTING "Build testsuite" ${PROJECT_IS_TOP_LEVEL})
option(${otr}_INSTALL_TESTING_LIBRARY "Install the compiled library used for testsuite" OFF)
message(STATUS "Showing option ${otr}_BUILD_TESTING: ${${otr}_BUILD_TESTING}")
message(STATUS "Showing option ${otr}_INSTALL_TESTING_LIBRARY: ${${otr}_INSTALL_TESTING_LIBRARY}")

# default to static libraries
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
    set(LIB_SUFFIX "64")

    # check whether 64-bit BLAS and LAPACK routines are called differently
    set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    check_fortran_function_exists("sgemm_64" BLAS_64)
    check_fortran_function_exists("sgesv_64" LAPACK_64)
    if(BLAS_64)
        add_compile_definitions("ddot=ddot_64")
        add_compile_definitions("dnrm2=dnrm2_64")
        add_compile_definitions("dgemv=dgemv_64")
        add_compile_definitions("dgemm=dgemm_64")
    endif()
    if(LAPACK_64)
        add_compile_definitions("dsyev=dsyev_64")
        add_compile_definitions("dsysv=dsysv_64")
        add_compile_definitions("zheev=zheev_64")
    endif()
else()
    set(LIB_SUFFIX "32")
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
    OUTPUT_NAME "opentrustregion_${LIB_SUFFIX}"
    EXPORT_NAME "opentrustregion${LIB_SUFFIX}"
    # library version
    VERSION ${PROJECT_VERSION}
    # interface version
    SOVERSION ${PROJECT_VERSION_MAJOR}
    # position-independent code - not automatically enabled for static libraries, but
    # required when linking a static library into a shared library
    POSITION_INDEPENDENT_CODE ON
)

# allow undefined symbols at link time on macOS for shared libraries
if(${BUILD_SHARED_LIBS} AND APPLE)
    set_target_properties(opentrustregion PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

# create target-specific version property
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
    add_library(otrtestsuite SHARED ${OPENTRUSTREGION_TESTS})

    # enable preprocessing
    set_source_files_properties(${OPENTRUSTREGION_TESTS}
                                PROPERTIES Fortran_PREPROCESS ON
    )

    # link against BLAS and LAPACK and against OpenTrustRegion
    target_link_libraries(otrtestsuite
                          PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
    )

    # install testsuite
    if (${otr}_INSTALL_TESTING_LIBRARY)
        install(TARGETS otrtestsuite
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${test_component} # Windows (.dll)
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${test_component} # shared libraries (.so, .dylib)
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${test_component} # static libraries (.a)
        )
    endif()
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
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${runtime_component} # Windows (.dll)
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${runtime_component} # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${runtime_component} # static libraries (.a)
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT ${runtime_component}
)

# install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT ${runtime_component}
)

# install targets
install(EXPORT library_interface
        FILE ${otr}Targets-${LIB_SUFFIX}.cmake
        NAMESPACE "${otr}::"
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
        COMPONENT ${runtime_component}
)

# install config files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}ConfigVersion.cmake
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
        COMPONENT ${runtime_component}
)
