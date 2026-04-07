cmake_minimum_required(VERSION 3.0)

project(cquadpack C)

include_directories(NumbaQuadpack/include)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(source_files
    NumbaQuadpack/src/dqag.c
    NumbaQuadpack/src/dqage.c
    NumbaQuadpack/src/dqagi.c
    NumbaQuadpack/src/dqagp.c
    NumbaQuadpack/src/dqags.c
    NumbaQuadpack/src/dqawc.c
    NumbaQuadpack/src/dqawce.c
    NumbaQuadpack/src/dqawf.c
    NumbaQuadpack/src/dqawfe.c
    NumbaQuadpack/src/dqawo.c
    NumbaQuadpack/src/dqaws.c
    NumbaQuadpack/src/dqawse.c
    NumbaQuadpack/src/dqc25c.c
    NumbaQuadpack/src/dqc25o.c
    NumbaQuadpack/src/dqc25s.c
    NumbaQuadpack/src/dqcheb.c
    NumbaQuadpack/src/dqext.c
    NumbaQuadpack/src/dqfour.c
    NumbaQuadpack/src/dqk15.c
    NumbaQuadpack/src/dqk15i.c
    NumbaQuadpack/src/dqk15w.c
    NumbaQuadpack/src/dqk21.c
    NumbaQuadpack/src/dqk31.c
    NumbaQuadpack/src/dqk41.c
    NumbaQuadpack/src/dqk51.c
    NumbaQuadpack/src/dqk61.c
    NumbaQuadpack/src/dqmomo.c
    NumbaQuadpack/src/dqng.c
    NumbaQuadpack/src/dqsort.c
    NumbaQuadpack/src/dqwgt.c
)
set(header_files
    NumbaQuadpack/include/cquadpack.h
)
set(BUILD_SHARED_LIBS ON)
add_library(cquadpack ${source_files} ${header_files})
message(STATUS "Will compile the library in shared mode")
set_property(TARGET cquadpack PROPERTY POSITION_INDEPENDENT_CODE 1)

# Add export header
include(GenerateExportHeader)
generate_export_header(cquadpack)
include_directories(${PROJECT_BINARY_DIR})

# Add the math libraries (gcc requires this)
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(cquadpack m)
endif()

set_target_properties(cquadpack PROPERTIES SUFFIX ".so")
install(TARGETS cquadpack DESTINATION ${CMAKE_SOURCE_DIR}/NumbaQuadpack/)

