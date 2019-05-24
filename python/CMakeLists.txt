cmake_minimum_required( VERSION 3.4.3 )

project( python CXX )

# ------------ Set up library for C++ bindings -------------------

# Directories to include header files from
include_directories( ../splinekernel/inc )

# Add shared library of python bindings
pybind11_add_module(pysplinekernel bindings/splinekernel_python.cpp )

target_link_libraries( pysplinekernel PRIVATE splinekernel )

# specify the relative path the shared library object shall be installed to
install( TARGETS pysplinekernel LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
  
# -------------- Set up copying python scipts ---------------------
set( PYTHON_SCRIPT_FOLDER scripts )

# find all python files in script folder
file( GLOB PYTHON_SCRIPTS ${PYTHON_SCRIPT_FOLDER}/*.py )

# add all files to install target
foreach( filename ${PYTHON_SCRIPTS} )
  install(FILES ${filename} DESTINATION ${CMAKE_INSTALL_PREFIX} )
endforeach()

