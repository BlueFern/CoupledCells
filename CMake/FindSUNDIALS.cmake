# - Try to find SUNDIALS
#

find_path (SUNDIALS_DIR include/sundials/sundials_config.h HINTS ENV SUNDIALS_DIR PATHS $ENV{HOME}/sundials DOC "Sundials Directory")

set(SUNDIALS_LIBRARIES)

IF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND YES)
  SET(SUNDIALS_INCLUDES ${SUNDIALS_DIR})
  find_path (SUNDIALS_INCLUDE_DIR sundials_config.h HINTS "${SUNDIALS_DIR}" PATH_SUFFIXES include/sundials NO_DEFAULT_PATH)
  list(APPEND SUNDIALS_INCLUDES ${SUNDIALS_INCLUDE_DIR})

 find_library(sundials_arkode NAMES        sundials_arkode HINTS "${SUNDIALS_DIR}/lib")
 find_library(sundials_nvecserial NAMES    sundials_nvecserial HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_cvode NAMES         sundials_cvode HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_cvodes NAMES        sundials_cvodes HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_fcvode NAMES        sundials_fcvode HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_fida NAMES          sundials_fida HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_fkinsol NAMES       sundials_fkinsol HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_fnvecparallel NAMES sundials_fnvecparallel HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_fnvecserial NAMES   sundials_fnvecserial HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_ida NAMES           sundials_ida HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_idas NAMES          sundials_idas HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_kinsol NAMES        sundials_kinsol HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_nvecparallel NAMES  sundials_nvecparallel HINTS "${SUNDIALS_DIR}/lib")
 #find_library(sundials_nvecserial NAMES    sundials_nvecserial HINTS "${SUNDIALS_DIR}/lib")
  
 set(SUNDIALS_LIBRARIES
 ${sundials_arkode}
 ${sundials_nvecserial}
 #${sundials_cvode}
 #${sundials_cvodes}
 #${sundials_fcvode}
 #${sundials_fida}
 #${sundials_fkinsol}
 #${sundials_fnvecparallel}
 #${sundials_fnvecserial}
 #${sundials_ida}
 #${sundials_idas}
 #${sundials_kinsol}
 #${sundials_nvecparallel}
 )
 	
ELSE(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)
  SET(SUNDIALS_FOUND NO)
  message(FATAL_ERROR "Cannot find SUNDIALS!")
ENDIF(EXISTS ${SUNDIALS_DIR}/include/sundials/sundials_config.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS DEFAULT_MSG SUNDIALS_LIBRARIES SUNDIALS_INCLUDES)
