#--------
#-search
#--------
include_directories(${PROJECT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR} )
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set( hepmc3search_sources
${PROJECT_SOURCE_DIR}/search/src/Relatives.cc
${PROJECT_SOURCE_DIR}/search/src/Selector.cc)

add_library(HepMC3search SHARED ${hepmc3search_sources})
target_link_libraries(HepMC3search  HepMC3)
target_include_directories(HepMC3search PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(HepMC3search  HepMC3)
if(${CMAKE_VERSION} VERSION_LESS "3.0.0")
  set_target_properties(HepMC3search PROPERTIES OUTPUT_NAME HepMC3search   SOVERSION 4
                                               ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/)

else()
  set_target_properties(HepMC3search PROPERTIES OUTPUT_NAME HepMC3search   SOVERSION 4
                                               ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)
endif()
# installs
install(TARGETS HepMC3search DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT searchlibs)
if (HEPMC3_BUILD_STATIC_LIBS)
  add_library(HepMC3search_static STATIC  ${hepmc3search_sources})
  target_include_directories(HepMC3search_static PUBLIC ${PROJECT_SOURCE_DIR}/include)
  target_compile_definitions(HepMC3search_static PRIVATE NO_DECLSPEC_StandardSelector)
  if(${CMAKE_VERSION} VERSION_LESS "3.0.0")
    set_target_properties(HepMC3search_static PROPERTIES OUTPUT_NAME HepMC3search-static
                                               ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/)
  else()
    set_target_properties(HepMC3search_static PROPERTIES OUTPUT_NAME HepMC3search-static
                                               ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
                                               LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>
                                               RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/${CMAKE_INSTALL_LIBDIR}/$<0:>)
  endif()
  install(TARGETS HepMC3search_static DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT searchlibs)
endif()

install(DIRECTORY ${PROJECT_BINARY_DIR}/search/ DESTINATION  ${CMAKE_INSTALL_LIBDIR}  COMPONENT searchlibs  FILES_MATCHING PATTERN "CMakeFiles"  EXCLUDE  PATTERN "src"  EXCLUDE PATTERN "include"  EXCLUDE)
install(DIRECTORY include/HepMC3 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT searchdevel)
