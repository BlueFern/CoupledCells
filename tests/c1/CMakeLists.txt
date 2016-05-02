
file(COPY ${CMAKE_SOURCE_DIR}/tests/c1/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c1/files DESTINATION "./")
add_test(NAME c1
	     COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1"
	     "${CMAKE_BINARY_DIR}/coupledCellsModel" 
	     -args 
	     -f config.txt 
	     -S solution
	     -T profiling
	     -t 10.0 
	     -w 1.0 
	     -i 1.e-2)
