project(espresso LANGUAGES C)

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB FILENAMES *.c)

# Surpress warnings from this external library
if (UNIX)
    set(CMAKE_C_FLAGS "-w")
endif()

add_library(espresso ${FILENAMES})
target_include_directories(espresso INTERFACE ${PROJECT_SOURCE_DIR})
