cmake_minimum_required(VERSION 3.16)
project(petsc_min LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Help FindPETSc locate your build (or pass these on the cmake command line)
if(NOT DEFINED ENV{PETSC_DIR})
  set(ENV{PETSC_DIR} "/home/arash/Projects/petsc")
endif()
if(NOT DEFINED ENV{PETSC_ARCH})
  set(ENV{PETSC_ARCH} "arch-linux-cuda")
endif()

find_package(PETSc REQUIRED)          # provides PETSC_INCLUDE_DIRS, PETSC_LIBRARIES, etc.

add_executable(petsc_hello main.cpp)
target_include_directories(petsc_hello PRIVATE ${PETSC_INCLUDE_DIRS})
target_link_libraries(petsc_hello PRIVATE ${PETSC_LIBRARIES})
# (Optional) warnings:
target_compile_options(petsc_hello PRIVATE -Wall -Wextra -Wpedantic)