configure_file(c256f_4x64_fitzroy.ll.in c256f_4x64_fitzroy.ll)
configure_file(c256f_8x32_fitzroy.ll.in c256f_8x32_fitzroy.ll)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c1/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c1/files DESTINATION "./")
add_custom_target(solution_c1 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c1 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
