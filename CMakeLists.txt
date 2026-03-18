cmake_minimum_required(VERSION 3.16)

project(simplex_method_solver VERSION 2.0 LANGUAGES CXX)

include(CheckIPOSupported)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(simplex_core
    src/simplex_solver.cpp
)

target_include_directories(simplex_core
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

if(MSVC)
    target_compile_options(simplex_core PRIVATE /W4 /permissive-)
else()
    target_compile_options(simplex_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(simplex
    src/simplex.cpp
)

target_link_libraries(simplex PRIVATE simplex_core)

if(MINGW)
    target_link_options(simplex PRIVATE -static -static-libgcc -static-libstdc++)
endif()

check_ipo_supported(RESULT SIMPLEX_IPO_SUPPORTED OUTPUT SIMPLEX_IPO_ERROR)
if(SIMPLEX_IPO_SUPPORTED)
    set_property(TARGET simplex_core simplex PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET simplex_core simplex PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
endif()

enable_testing()

add_executable(simplex_tests
    tests/simplex_tests.cpp
)

target_link_libraries(simplex_tests PRIVATE simplex_core)

if(MINGW)
    target_link_options(simplex_tests PRIVATE -static -static-libgcc -static-libstdc++)
endif()

add_test(NAME simplex_tests COMMAND simplex_tests)
