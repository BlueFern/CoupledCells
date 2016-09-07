configure_file(c32_4procs_fitzroy.ll.in c32_4procs_fitzroy.ll)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c32_4procs/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c32_4procs/files DESTINATION "./")
add_custom_target(solution_c32_4procs ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c32_4procs ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
