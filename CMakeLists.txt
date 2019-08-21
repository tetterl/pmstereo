cmake_minimum_required(VERSION 3.4)

project(pmstereo CXX)

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -mavx2 -mfma")

# Compile fast_exp just to get the assembly output for inspection. It's slower to actually use this library
# rather than inlining though
add_library(fast_exp_instr STATIC fast_exp_instr.cpp)
set_target_properties(fast_exp_instr PROPERTIES COMPILE_FLAGS "-save-temps")
add_library(fast_exp_array STATIC fast_exp_array.cpp)
set_target_properties(fast_exp_array PROPERTIES COMPILE_FLAGS "-save-temps")
add_library(fast_exp_array_decl STATIC fast_exp_array_declaration.cpp)
set_target_properties(fast_exp_array_decl PROPERTIES COMPILE_FLAGS "-save-temps")

add_executable(pm pm.cpp)
target_link_libraries(pm ${PNG_LIBRARIES})

# save asssembly file
set_target_properties(pm PROPERTIES COMPILE_FLAGS "-save-temps")

add_executable(eval eval.cpp)
target_link_libraries(eval ${PNG_LIBRARIES})
