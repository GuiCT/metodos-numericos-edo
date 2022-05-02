cmake_minimum_required (VERSION 3.8)

project (NumericalMethods)

option(USE_AVX "Build using AVX2 instructions" ON)

#[[Biblioteca:
Método numérico de CashKarp#]]

add_library(CashKarp STATIC
    ${PROJECT_SOURCE_DIR}/CashKarp/CashKarp.cpp
)

target_include_directories(CashKarp PUBLIC
    ${PROJECT_SOURCE_DIR}/CashKarp
)

#[[Arquivo de execução, testando métodos]]

add_executable(NumericalMethods 
    ${PROJECT_SOURCE_DIR}/main.cpp
)

if(USE_AVX)
    target_compile_options(CashKarp PUBLIC /arch:AVX2)
endif(USE_AVX)

#[[Adicionando biblioteca no projeto]]

target_link_libraries(NumericalMethods PRIVATE
    CashKarp
)