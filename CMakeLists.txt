cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

project(CoupledCells)

set(PROJECT_NAME "Coupled Arterial Cells")
set(PROJECT_COPYRIGHT "UC HPC, University of Canterbury, NZ")

# TODO: Update version numbers.
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 9)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(PROJECT_DESCRIPTION "Large-scale simulation of chemical species transport in coupled arterial cells.")

# Add local CMake directory to CMake module path to load project-specific modules. 
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Take care of the non-standard location of malloc includes on Mac.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Unable to run without MPI.
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# Unable to run without HDF5.
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})

# ODE implementation choice.
set(ODE_SOLVER "RK_suite" CACHE STRING
	"ODE solver implementation to be used in the simulations.")

# Available options for ODE implementation.
set(ODE_SOLVER_CHOICES "RK_suite;SUNDIALS_arkode;BOOST_odeint" CACHE INTERNAL
	"List of values for choosing a solver to be compiled against.")

# Enable ODE choices in the GUI drop-down menue.
# The appropriate choice will have to be typed in when using CLI.  
set_property(CACHE ODE_SOLVER PROPERTY STRINGS ${ODE_SOLVER_CHOICES})

# In this case we want to compile RK Suite in the ext/rksuite-1.0 directory.
if(${ODE_SOLVER} STREQUAL "RK_suite")
	message(STATUS "Selected ODE solver: Local RK Suite")
	
	# Pass preprocessor macro definition to the compiler.
	add_definitions(-DRK_SUITE)
	
	# Add RK Suite library.
	include_directories(${PROJECT_SOURCE_DIR}/ext/rksuite-1.0)
	add_subdirectory(${PROJECT_SOURCE_DIR}/ext/rksuite-1.0)
	
	# Remember the library for the ODE solver implementation to link in to the executable.
	set(ODE_SOLVER_LIBS RKSuite)
	
	set(ODE_SOLVER_WRAPPER src/solver_rk.cpp)

# In this case we want to find the system SUNDIALS libraries, including ARKODE.
elseif(${ODE_SOLVER} STREQUAL "SUNDIALS_arkode")
	message(STATUS "Selected ODE solver: system SUNDIALS arkode")
	
	# Passe preprocessor macro definition to the compiler.
	add_definitions(-DARK_ODE)
	
	# Find SUNDIALS
	find_package(SUNDIALS REQUIRED)
	# Add the location of the SUNDIALS headers.
	include_directories(${SUNDIALS_INCLUDES})
	
	# Remember the libraries for the ODE solver implementation to link in to the executable.
	set(ODE_SOLVER_LIBS ${SUNDIALS_LIBRARIES})
	
	set(ODE_SOLVER_WRAPPER src/solver_arkode.cpp)
	
# In this case we want to find the system BOOST libraries, including ODEINT.
elseif(${ODE_SOLVER} STREQUAL "BOOST_odeint")
	message(STATUS "Selected ODE solver: system Boost odeint")
	
	# Passe preprocessor macro definition to the compiler.
	add_definitions(-DBOOST_ODEINT)
	
	# Find Boost
	find_package(Boost REQUIRED)
	# Add the location of the Boost headers.
	include_directories(${Boost_INCLUDE_DIRS})
	
	set(ODE_SOLVER_WRAPPER src/solver_odeint.cpp)

endif(${ODE_SOLVER} STREQUAL "RK_suite")

file(GLOB HDRs src/*.h)

set (SRCs 
	src/checkpoint.cpp 
	src/koenigsberger_model.cpp
	src/tsoukias_smc_multicell.cpp
	src/update.cpp
	src/compute.cpp
	src/main.cpp 
	src/writeHDF5.cpp
	src/debug_file.cpp 
	src/topology.cpp
	src/gather.cpp
	${ODE_SOLVER_WRAPPER})

add_executable(coupledCellsModel ${SRCs} ${HDRs})

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

target_link_libraries(coupledCellsModel 
	${MPI_LIBRARIES}
	${HDF5_LIBRARIES}
	${ODE_SOLVER_LIBS})

find_package(Doxygen)

find_Package(Graphviz)
if(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "YES")
	set(CALL_GRAPH "YES")
	set(CALLER_GRAPH "YES")
else(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "NO")
	set(CALL_GRAPH "NO")
	set(CALLER_GRAPH "NO")
endif(DOXYGEN_DOT_EXECUTABLE)

set(DOXY_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/doxygen)
file(MAKE_DIRECTORY ${DOXY_OUTPUT_DIR})

set(DOXY_API_LINK_TAGS ${DOXY_OUTPUT_DIR}/DoxyLink.tag)

configure_file(
	${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.in
	${PROJECT_BINARY_DIR}/Doxyfile
	@ONLY)

if(DOXYGEN_FOUND)
	add_custom_target(doxyDocs
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

find_package(Sphinx)

if(NOT DEFINED SPHINX_THEME)
	set(SPHINX_THEME default)
endif()
if(NOT DEFINED SPHINX_THEME_DIR)
	set(SPHINX_THEME_DIR)
endif()

set(SPHINX_BUILD_DIR ${PROJECT_BINARY_DIR}/_sphinxbuild)
# Muffle the warning about missing directory.
# Might have to use _static directory in source to copy files here if necessary.
file(MAKE_DIRECTORY ${SPHINX_BUILD_DIR}/_static)
set(SPHINX_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/sphinx)

configure_file(
    ${PROJECT_SOURCE_DIR}/doc/sphinx/conf.py.in
    ${SPHINX_BUILD_DIR}/conf.py
    @ONLY)

if(SPHINX_FOUND)
	add_custom_target(sphinxDocs
    COMMAND ${SPHINX_EXECUTABLE}
        -q -b html
        -c ${SPHINX_BUILD_DIR}
        -d ${SPHINX_BUILD_DIR}
        ${PROJECT_SOURCE_DIR}/doc/sphinx
        ${SPHINX_OUTPUT_DIR}
    COMMENT "Building HTML documentation with Sphinx" VERBATIM)
endif(SPHINX_FOUND)

