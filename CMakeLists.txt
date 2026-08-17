cmake_minimum_required(VERSION 3.20)
project(deformation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only set a default build type for single-config generators (Makefiles/Ninja).
# Multi-config generators (e.g. Visual Studio) ignore this and use --config instead.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(OpenMP)

# Eigen
include(FetchContent)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(eigen)

# polyscope vendors its own glfw/glad, so no separate windowing dependency is needed
add_subdirectory(deps/polyscope)

file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

list(REMOVE_ITEM PROJECT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

add_library(cd_core ${PROJECT_SOURCES})
target_include_directories(cd_core PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(cd_core PUBLIC
  Eigen3::Eigen
  polyscope
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(cd_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(cd_core PUBLIC PM_USE_OPENMP=1)
else()
  message(WARNING "OpenMP not found. Building without OpenMP parallelism.")
endif()

add_executable(profile_mover src/main.cpp)
target_link_libraries(profile_mover PRIVATE cd_core)