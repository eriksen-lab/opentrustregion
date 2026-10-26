# builds a C and a Fortran program against an installed OpenTrustRegion, which the CI
# job of the CMake installation runs, so that the header, the Fortran modules and the
# libraries the package links are checked to reach programs that use it
cmake_minimum_required(VERSION 3.24)
project(consumer LANGUAGES C Fortran)

# integer size of the installed library in bits
set(INTEGER_BITS 32 CACHE STRING "Integer size of the installed library in bits")

find_package(OpenTrustRegion CONFIG REQUIRED COMPONENTS int${INTEGER_BITS})

add_executable(consume_c consume.c)
target_link_libraries(consume_c PRIVATE OpenTrustRegion::OpenTrustRegion)

add_executable(consume_fortran consume.f90)
target_link_libraries(consume_fortran PRIVATE OpenTrustRegion::OpenTrustRegion)

# print the properties of the imported targets for the CI log
include(CMakePrintHelpers)
cmake_print_properties(
  TARGETS
    OpenTrustRegion::opentrustregion${INTEGER_BITS}
    OpenTrustRegion::opentrustregion
    OpenTrustRegion::OpenTrustRegion
  PROPERTIES
    IMPORTED_LINK_DEPENDENT_LIBRARIES
    IMPORTED_LOCATION
    INTERFACE_LOCATION
    INTERFACE_LINK_LIBRARIES
    INTERFACE_INCLUDE_DIRECTORIES
    INTERFACE_COMPILE_DEFINITIONS
    INTERFACE_COMPILE_OPTIONS
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
  )
