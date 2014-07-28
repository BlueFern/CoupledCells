# - Try to find SUNDIALS
#

find_path (SUNDIALS_DIR include/sundials/sundials_config.h HINTS ENV SUNDIALS_DIR PATHS $ENV{HOME}/sundials DOC "Sundials Directory")

IF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND YES)
  SET(SUNDIALS_INCLUDES ${SUNDIALS_DIR})
  find_path (SUNDIALS_INCLUDE_DIR sundials_config.h HINTS "${SUNDIALS_DIR}" PATH_SUFFIXES include/sundials NO_DEFAULT_PATH)
  list(APPEND SUNDIALS_INCLUDES ${SUNDIALS_INCLUDE_DIR})
  find_library(SUNDIALS_LIBRARIES NAMES
	sundials_cvodes
	sundials_fcvode
	sundials_fida
	sundials_fkinsol
	sundials_fnvecparallel
	sundials_fnvecserial
	sundials_ida
	sundials_idas
	sundials_kinsol
	sundials_nvecparallel
	sundials_nvecserial
 	HINTS "${SUNDIALS_DIR}/lib")
ELSE(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND NO)
  message(FATAL_ERROR "Cannot find SUNDIALS!")
ENDIF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS DEFAULT_MSG SUNDIALS_LIBRARIES SUNDIALS_INCLUDES)