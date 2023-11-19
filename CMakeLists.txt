cmake_minimum_required(VERSION 3.10)
project(QuadRotorSim VERSION 0.0.0)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER -O3")

# Specify the source files
file(GLOB_RECURSE SOURCES src/*.cpp)

add_subdirectory(3rdparty/matplotplusplus)

# Find Armadillo, BLAS, LAPACK, Boost and OpenMP
find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)



# Include directories for headers
include_directories(
    include
    ${ARMADILLO_INCLUDE_DIRS}
    ${BLAS_INCLUDE_DIRS}
    ${LAPACK_INCLUDE_DIRS}
    ${OpenMP_INCLUDE_DIRS}
)

# Configure version.h
configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.in
    ${CMAKE_SOURCE_DIR}/include/version.h
)

# Generate version.h
message(STATUS "QuadRotorSim version: ${PROJECT_VERSION}")

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link Armadillo library
target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} OpenMP::OpenMP_CXX matplot)

# Set the output directory for binary files
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Set the output directory for libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)