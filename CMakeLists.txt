cmake_minimum_required (VERSION 2.6)

project(CoupledCells)

set(PROJECT_NAME "Coupled Arterial Cells")
set(PROJECT_COPYRIGHT "BlueFern, University of Canterbury, NZ")

set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 9)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(PROJECT_DESCRIPTION "Large-scale simulation of chemical species transport in coupled arterial cells.")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

find_package(SUNDIALS REQUIRED)
include_directories(${SUNDIALS_INCLUDES})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})

file(GLOB HDRs src/*.h)
file(GLOB SRCs src/*.cpp)

include_directories(${PROJECT_SOURCE_DIR}/contrib)
add_subdirectory(contrib)

add_executable(coupledCellsModel ${SRCs} ${HDRs})

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

target_link_libraries(coupledCellsModel ${SUNDIALS_LIBRARIES} ${HDF5_LIBRARIES} ${MPI_LIBRARIES} RKSuite)

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
# Muffle a warning about missing directory.
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

