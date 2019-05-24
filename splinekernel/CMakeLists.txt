cmake_minimum_required( VERSION 3.4.3 )

project( splinekernel CXX )

# Directories to include header files from
include_directories( inc )

# Gather source files into CMake variables
file( GLOB SOURCE_FILES src/*.cpp)
file( GLOB TEST_FILES test/*.cpp)
file( GLOB HEADER_FILES inc/*.h* )

# This enables exporting all symbols to the dll on windows
set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )

# ----------------- Set up spline library -------------------------

# Add library from the collected source files. The headers are given so visual studio displays them
add_library( splinekernel SHARED ${SOURCE_FILES} ${HEADER_FILES} ) 

# specify the relative path the shared library object shall be installed to
if( WIN32 )
  install( TARGETS splinekernel RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX} )
else( )
  install( TARGETS splinekernel LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
endif( )

# ------------------- Set up unit tests ---------------------------

# Add an executable and link to the library
add_executable( splinekernel_testrunner ${TEST_FILES} )

target_link_libraries( splinekernel_testrunner splinekernel )

# specify the relative path the testrunner shall be installed to
install( TARGETS splinekernel_testrunner RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX} )
