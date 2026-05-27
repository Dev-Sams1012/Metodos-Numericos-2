cmake_minimum_required(VERSION 3.15)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
project(NumericalMethods VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --------------------------------------------------------------------------
# Dependência externa: muParser (parser de expressões matemáticas)
# --------------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
    muparser
    GIT_REPOSITORY https://github.com/beltoforion/muparser.git
    GIT_TAG        v2.3.4
)

# === A SOLUÇÃO ESTÁ AQUI ===
# Força a compilação como biblioteca estática para evitar erro de DLL no Windows
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
# ===========================

set(ENABLE_SAMPLES  OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP   OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(muparser)

# --------------------------------------------------------------------------
# Coleta de fontes (src/ e utils/)
# --------------------------------------------------------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/utils/*.cpp"
)

# --------------------------------------------------------------------------
# Executável principal
# --------------------------------------------------------------------------
add_executable(numerical_methods
    main.cpp
    ${SOURCES}
)

target_link_options(numerical_methods PRIVATE -mconsole)

target_include_directories(numerical_methods PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/utils
)

target_link_libraries(numerical_methods PRIVATE muparser)

# Flags de compilação
target_compile_options(numerical_methods PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2 -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /W4>
)

message(STATUS "NumericalMethods configurado com sucesso!")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")