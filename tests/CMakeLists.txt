add_test(NAME c1
	     COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1"
	     "${CMAKE_BINARY_DIR}/coupledCellsModel" 
	     -args "-f ${CMAKE_SOURCE_DIR}/tests/config_c1.txt -S solution_c1 -T profiling_c1 -t 10.0 -w 1.0 -i 1.e-2")
