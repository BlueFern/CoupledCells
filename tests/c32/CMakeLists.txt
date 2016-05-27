configure_file(c32_fitzroy.ll.in c32_fitzroy.ll)
configure_file(c32_compare_results.sh.in c32_compare_results.sh)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c32/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c32/files DESTINATION "./")
add_custom_target(solution_c32 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c32 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
add_test(NAME c32
	     COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "32"
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
    add_test(NAME c1_compare_results
             COMMAND "sh" "c1_compare_results.sh")
    set_tests_properties(c1_compare_results
             PROPERTIES FAIL_REGULAR_EXPRESSION 
             "[1-9][0-9]*[ ]*differences found")
endif()
