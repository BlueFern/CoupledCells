configure_file(c1008_16x64_fitzroy.ll.in c1008_16x64_fitzroy.ll)
configure_file(c1008_bgp.ll.in c1008_bgp.ll)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c1008/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c1008/files DESTINATION "./")
add_custom_target(solution_c1008 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c1008 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
