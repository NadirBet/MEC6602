cmake_minimum_required(VERSION 3.16)
project(EulerOGRIDSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (MSVC)
    add_compile_options(/MP /O2 /permissive-)
    add_compile_definitions(_USE_MATH_DEFINES)
else()
    add_compile_options(-O3 -Wall -Wextra -Wpedantic)
endif()

add_executable(euler_solver
    Solver.cpp
)



# So both can find Mesh.hpp, Flux.hpp, Boundary.hpp, Types.hpp, etc.
target_include_directories(euler_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
