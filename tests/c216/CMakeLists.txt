configure_file(c216_fitzroy.ll.in c216_fitzroy.ll)
configure_file(c216_compare_results.sh.in c216_compare_results.sh)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c216/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c216/files DESTINATION "./")
add_custom_target(solution_c216 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c216 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
add_test(NAME c216
	     COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "216"
	     "${CMAKE_BINARY_DIR}/coupledCellsModel" 
	     -f config.txt 
	     -S solution
	     -T profiling
	     -t 100.0 
	     -w 1.0 
	     -i 1.e-2
             -C 1
	     -R 0)

# success is defined as the absence of the regular expression
# "[1-9][0-9]*[ ]*differences found" in the output
if ("${HAVE_REF_RESULTS}" STREQUAL "YES")
    add_test(NAME c32_compare_results
             COMMAND "sh" "c216_compare_results.sh")
    set_tests_properties(c216_compare_results
             PROPERTIES FAIL_REGULAR_EXPRESSION 
             "[1-9][0-9]*[ ]*differences found")
endif()
