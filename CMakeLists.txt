cmake_minimum_required(VERSION 3.20)
project(GodunovSolver)  # ← Измените название проекта

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI QUIET COMPONENTS CXX)

# Добавьте все исходные файлы
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/parser.cpp
    src/point.cpp
    src/grid.cpp
    src/solver.cpp
    src/flic.cpp
    src/mader.cpp
    src/script.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
)
if(MSVC)
    add_compile_options(/O2 /fp:fast /arch:AVX2)
else()
    add_compile_options(-O3 -ffast-math -march=native)
endif()

if(MPI_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
endif()

# Уберите эту строку - она не нужна, если все файлы уже добавлены выше
# add_subdirectory("src")
