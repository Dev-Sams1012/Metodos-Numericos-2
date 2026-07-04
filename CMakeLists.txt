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

# Remove o proc_imgs.cpp da lista global para não dar conflito de 'main' duplicado
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/proc_imgs.cpp")

# --------------------------------------------------------------------------
# Executável 1: Métodos Numéricos Principal (CLI)
# --------------------------------------------------------------------------
add_executable(numerical_methods
    main.cpp
    ${SOURCES}
)

if(MINGW)
    target_link_options(numerical_methods PRIVATE -mconsole)
endif()

target_include_directories(numerical_methods PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/utils
)

target_link_libraries(numerical_methods PRIVATE muparser)

# --------------------------------------------------------------------------
# Executável 2: Processador de Imagens (Detecção de Bordas)
# --------------------------------------------------------------------------
add_executable(image_processor
    atividades/proc_imgs.cpp
    ${SOURCES} # Inclui as classes do projeto (Matrix.cpp, FDMethod1D.cpp, etc.)
)

if(MINGW)
    target_link_options(image_processor PRIVATE -mconsole)
endif()

# Garante que o image_processor veja as pastas locais e os arquivos gerados do muParser
target_include_directories(image_processor PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/utils
    ${muparser_SOURCE_DIR}/include  # <-- Força o caminho do muParser baixado
)

# Linka a biblioteca do muParser para puxar as dependências e headers corretamente
target_link_libraries(image_processor PRIVATE muparser)

# Aplica as mesmas flags de compilação otimizadas
target_compile_options(image_processor PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2 -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /W4>
)

target_compile_options(numerical_methods PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2 -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /W4>
)

message(STATUS "NumericalMethods configurado com sucesso!")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")