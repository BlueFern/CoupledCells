configure_file(c256f_4x64_fitzroy.ll.in c256f_4x64_fitzroy.ll)
configure_file(c256f_8x32_fitzroy.ll.in c256f_8x32_fitzroy.ll)
configure_file(c256f_bgp.ll.in c256f_bgp.ll)
configure_file(c256f_compare_results.sh.in c256f_compare_results.sh)
configure_file(c256f_create_slog2.sh.in c256f_create_slog2.sh)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c256f/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c256f/files DESTINATION "./")
add_custom_target(solution_c256f ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c256f ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")

if ("${HAVE_REF_RESULTS}" STREQUAL "YES")
    add_test(NAME c256f_compare_results
             COMMAND "sh" "c256f_compare_results.sh")
    set_tests_properties(c256f_compare_results
             PROPERTIES FAIL_REGULAR_EXPRESSION 
             "[0-9][ ]*differences found")
endif()
