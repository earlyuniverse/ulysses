
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was HepMC3Config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

SET(HEPMC3_VERSION 3.02.04)
SET(HEPMC3_VERSION_MAJOR  3)
SET(HEPMC3_VERSION_MINOR  2)
SET(HEPMC3_VERSION_PATCH  4)


set_and_check(HEPMC3_INCLUDE_DIR ${PACKAGE_PREFIX_DIR}/include)
set(HEPMC3_CXX_STANDARD 11)

if(EXISTS ${PACKAGE_PREFIX_DIR}/share/HepMC3/interfaces)
  set(HEPMC3_INTERFACES_DIR ${PACKAGE_PREFIX_DIR}/share/HepMC3/interfaces)
endif()

find_library(HEPMC3_LIB NAMES HepMC3 HINTS ${PACKAGE_PREFIX_DIR}/lib)
find_library(HEPMC3_SEARCH_LIB NAMES HepMC3search HINTS ${PACKAGE_PREFIX_DIR}/lib)
find_library(HEPMC3_ROOTIO_LIB NAMES HepMC3rootIO HINTS ${PACKAGE_PREFIX_DIR}/lib)

set(HEPMC3_LIBRARIES ${HEPMC3_LIB})
if(EXISTS ${HEPMC3_SEARCH_LIB})
  list( APPEND  HEPMC3_LIBRARIES ${HEPMC3_SEARCH_LIB})
endif()
if(EXISTS ${HEPMC3_ROOTIO_LIB})
  list( APPEND  HEPMC3_LIBRARIES ${HEPMC3_ROOTIO_LIB})
endif()
