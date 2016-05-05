cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

project(CoupledCells)

set(PROJECT_NAME "Coupled Arterial Cells")
set(PROJECT_COPYRIGHT "UC HPC, University of Canterbury, NZ")

set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(PROJECT_DESCRIPTION "Large-scale simulation of chemical species transport in coupled arterial cells.")

# Build type defaults to Release unless specified
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif (NOT CMAKE_BUILD_TYPE)

# Add local CMake directory to CMake module path to load project-specific modules. 
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Take care of the non-standard location of malloc includes on Mac.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Unable to run without MPI.
if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        # AIX (fitzroy) needs some help
        set(POE "/usr/lpp/ppe.poe")
        if(EXISTS "${POE}/include" AND EXISTS "${POE}/lib")
                set(MPI_INCLUDE_PATH "${POE}/include")
                set(MPI_LIBRARIES "${POE}/lib/libmpi_r.a;${POE}/lib/libvtd_r.a")
                set(MPI_C_INCLUDE_PATH "${MPI_INCLUDE_PATH}")
                set(MPI_C_LIBRARIES "${MPI_LIBRARIES}")
                set(MPI_CXX_INCLUDE_PATH "${MPI_INCLUDE_PATH}")
                set(MPI_CXX_LIBRARIES "${MPI_LIBRARIES}")
                message(STATUS "Setting MPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}")
                message(STATUS "Setting MPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}")
		if (CMAKE_BUILD_TYPE STREQUAL "Release")
		        set(ADD_FLAGS "-q64 -O3 -qstrict -qmaxmem=-1 -qtune=pwr6 -qarch=pwr6 -qcache=auto -qhot -qsimd -qipa")
		        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_FLAGS}")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_FLAGS}")
			set(CMAKE_Fortran_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_FLAGS}")
		endif (CMAKE_BUILD_TYPE STREQUAL "Release")
        endif(EXISTS "${POE}/include" AND EXISTS "${POE}/lib")
endif(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

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
	
	# Pass preprocessor macro definition to the compiler.
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
	
	# Pass preprocessor macro definition to the compiler.
	add_definitions(-DBOOST_ODEINT)
	
	# Find Boost
	find_package(Boost REQUIRED)
	# Add the location of the Boost headers.
	include_directories(${Boost_INCLUDE_DIRS})
	
	set(ODE_SOLVER_WRAPPER src/solver_odeint.cpp)

endif(${ODE_SOLVER} STREQUAL "RK_suite")

# Optionally turn on performance monitoring with TAU on AIX
# if ${TAU_MAKEFILE} is set
if (${CMAKE_SYSTEM_NAME} STREQUAL "AIX" AND (NOT "${TAU_MAKEFILE}" STREQUAL ""))
    get_filename_component(TAU_LIBRARY_DIR "${TAU_MAKEFILE}" DIRECTORY)
    get_filename_component(TAU_DIR "${TAU_LIBRARY_DIR}" DIRECTORY)
    set(ENV{TAU_MAKEFILE} "${TAU_MAKEFILE}")
    set(ENV{PATH} "${TAU_DIR}/bin:$ENV{PATH}")
    message(STATUS "**********************************************************")
    #message("PATH: $ENV{PATH}")
    message(STATUS "Please set the follwing prior to typing 'make':")
    message(STATUS "export PATH=${TAU_DIR}/bin:$PATH")
    message(STATUS "export TAU_MAKEFILE=$ENV{TAU_MAKEFILE}")
    message(STATUS "**********************************************************")
    set(CMAKE_C_COMPILER "tau_cc.sh")
    set(CMAKE_C_FLAGS "-optVerbose -O2 -qstrict -qmaxmem=-1 -qarch=pwr6 -qtune=pwr6")
    set(CMAKE_CXX_COMPILER "tau_cxx.sh")
    set(CMAKE_CXX_FLAGS "-optVerbose -O2 -qstrict -qmaxmem=-1 -qarch=pwr6 -qtune=pwr6")
    set(CMAKE_Fortran_COMPILER "tau_f90.sh") 
    set(CMAKE_Fortran_FLAGS "-optVerbose ${CMAKE_Fortran_FLAGS} -qstrict -qfixed -O2 -qmaxmem=-1 -qarch=pwr6 -qtune=pwr6")
endif()

file(GLOB HDRs src/*.h)

set (SRCs 
	src/time_profiling.cpp 
	src/koenigsberger_model.cpp
	src/tsoukias_model.cpp
	src/update.cpp
	src/compute.cpp
	src/simulation_main.cpp 
	src/writeHDF5.cpp
	src/debug_helpers.cpp 
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

message(STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS}")

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

# CTest
enable_testing()

add_subdirectory(tests)