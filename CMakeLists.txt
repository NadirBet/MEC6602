cmake_minimum_required(VERSION 3.16)
project(EulerOGRIDSolver LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Path to Eigen
include_directories("C:/Users/nadir/Desktop/CPP/eigen-3.4.0")

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Collect sources for main solver ---
set(SOURCES
    main.cpp
    FASMultigrid.cpp          # <-- new multigrid implementation
)

set(HEADERS
    mesh.hpp
    Flux.hpp
    Boundary.hpp
    Types.hpp
    ImplicitResidualSmoothing.hpp
    FASMultigrid.hpp          # <-- multigrid header
)

# -- Main solver executable --
add_executable(euler_solver ${SOURCES} ${HEADERS})

target_include_directories(euler_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    target_compile_options(euler_solver PRIVATE /W4 /O2 /permissive- /MP)
    target_compile_definitions(euler_solver PRIVATE _USE_MATH_DEFINES)
else()
    target_compile_options(euler_solver PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

# -- Cartesian constant-state wall test executable --
# Uses TEST_CartesianWall.cpp as its own main()
add_executable(cartesian_wall_test
    TEST_CartesianWall.cpp
    ${HEADERS}
)

target_include_directories(cartesian_wall_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    target_compile_options(cartesian_wall_test PRIVATE /W4 /O2 /permissive- /MP)
    target_compile_definitions(cartesian_wall_test PRIVATE _USE_MATH_DEFINES)
else()
    target_compile_options(cartesian_wall_test PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

# --- Test runner (optional) ---
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_runner.cpp)
    add_executable(tests_runner
        test_runner.cpp
        ${HEADERS}
    )
    target_include_directories(tests_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    if(MSVC)
        target_compile_options(tests_runner PRIVATE /W4 /O2 /permissive- /MP)
    else()
        target_compile_options(tests_runner PRIVATE -O3 -Wall -Wextra -Wpedantic)
    endif()
endif()
