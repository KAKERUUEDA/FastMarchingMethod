set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set the C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.0)

# Set the compiler to GCC
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})

project(${ProjectId})

file(GLOB project_sources src/*.cpp)
file(GLOB project_headers lib/*.h)

add_executable(${ProjectId} ${project_sources} ${project_headers})
