cmake_minimum_required(VERSION 3.17)
project(ARC_Simulation_Reconstruction
	VERSION 1.0
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT 6.22 CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

add_compile_options(-Wall)
add_compile_options(-Wunused)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-O3)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

add_subdirectory(${CMAKE_SOURCE_DIR}/apps)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

target_compile_definitions(ARC_Simulation_Reconstruction PUBLIC INSTALL_DIR="${CMAKE_CURRENT_SOURCE_DIR}/build/")
