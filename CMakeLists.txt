cmake_minimum_required (VERSION 2.6)

project(CoupledCells)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

find_package(SUNDIALS REQUIRED)
include_directories(${SUNDIALS_INCLUDES})

file(GLOB HDRs src/*.h)
file(GLOB SRCs src/*.cpp)

add_executable(coupledCellsModel ${SRCs} ${HDRs})

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}" FORCE)

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}" FORCE)

target_link_libraries(coupledCellsModel ${SUNDIALS_LIBRARIES} ${MPI_LIBRARIES})

find_package(Doxygen)

set(DOXY_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc)

if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/doc)
	configure_file(
	${PROJECT_SOURCE_DIR}/Doxyfile.in
	${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doxyDoc
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)
