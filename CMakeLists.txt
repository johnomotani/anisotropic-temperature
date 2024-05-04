cmake_minimum_required(VERSION 3.18)
cmake_policy(VERSION 3.18)

project(anisotropic-temperature LANGUAGES CXX)
find_package(bout++ REQUIRED)
add_executable(anisotropic-temperature anisotropic-temperature.cxx)
target_link_libraries(anisotropic-temperature PRIVATE bout++::bout++)

# Copy data directory with input file to build directory
add_custom_command(TARGET anisotropic-temperature POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data $<TARGET_FILE_DIR:anisotropic-temperature>/data)

# Copy useful scripts, etc. into build directory
add_custom_command(TARGET anisotropic-temperature POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/*.py $<TARGET_FILE_DIR:anisotropic-temperature>)
