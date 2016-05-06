configure_file(c256f_4x64_fitzroy.ll.in c256f_4x64_fitzroy.ll)
configure_file(c256f_8x32_fitzroy.ll.in c256f_8x32_fitzroy.ll)
configure_file(c256f_bgp.ll.in c256f_bgp.ll)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c256f/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c256f/files DESTINATION "./")
add_custom_target(solution_c256f ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c256f ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
