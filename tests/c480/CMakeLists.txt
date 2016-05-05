configure_file(c480_4x64_fitzroy.ll.in c480_4x64_fitzroy.ll)
configure_file(c480_8x32_fitzroy.ll.in c480_8x32_fitzroy.ll)

file(COPY ${CMAKE_SOURCE_DIR}/tests/c480/config.txt DESTINATION "./")
file(COPY ${CMAKE_SOURCE_DIR}/tests/c480/files DESTINATION "./")
add_custom_target(solution_c480 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "solution")
add_custom_target(profiling_c480 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "profiling")
