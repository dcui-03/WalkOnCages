# cmake -S . -B build   -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++   -DOpenMP_ROOT=$(brew --prefix libomp)
cmake_minimum_required(VERSION 3.20)
project(deformation LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Only set default build type for single-config generators
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
  # set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


#if(APPLE)
#    set(CMAKE_OSX_ARCHITECTURES "arm64")
#endif()

# Libigl
set(LIBIGL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include)
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
# Other modules you could enable
#igl_include(embree)
#igl_include(imgui)
#igl_include(opengl)
#igl_include(stb)
#igl_include(predicates)
#igl_include(xml)
#igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
#igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

add_subdirectory(deps/polyscope)

file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

list(REMOVE_ITEM PROJECT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

add_library(sbc_core ${PROJECT_SOURCES})
target_include_directories(sbc_core PUBLIC
  #"${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

# Link igl (and the glfw module) to your project
target_link_libraries(sbc_core PUBLIC 
  igl::glfw
  ## Other modules you could link to
  # igl::embree
  # igl::imgui
  # igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  # igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  #igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  #igl_restricted::triangle
  polyscope
  #geometry-central
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(sbc_core PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(sbc_core PUBLIC PM_USE_OPENMP=1)
else()
  message(WARNING "OpenMP not found. Building without OpenMP parallelism.")
endif()

# Add your project files
add_executable(profile_mover src/main.cpp)

# 
target_link_libraries(profile_mover PRIVATE sbc_core)