cmake_minimum_required(VERSION 3.14)
project(gsea CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Eigen3_DIR /opt/homebrew/Cellar/eigen/5.0.0)

# Find Eigen
find_package(Eigen3 5.0 REQUIRED NO_MODULE)

# Source files
set(SOURCES
        src/main.cpp
        src/types/expression_data.cpp
        src/types/sample_data.cpp
        src/types/gene_set.cpp
        src/data_loader/expression_loader.cpp
        src/data_loader/sample_loader.cpp
        src/data_loader/geneset_loader.cpp
        src/gsea/ranking.cpp
        src/gsea/enrichment.cpp
        src/gsea/statistics.cpp
        src/gsea/analyzer.cpp
)

# Executable
add_executable(gsea ${SOURCES})

# Include directories
target_include_directories(gsea PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

# Link libraries
target_link_libraries(gsea Eigen3::Eigen)

# Compiler options
if(MSVC)
    target_compile_options(gsea PRIVATE /W4)
else()
    target_compile_options(gsea PRIVATE -Wall -Wextra -O2 -march=native)
endif()