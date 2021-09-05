cmake_minimum_required(VERSION 3.0)

project(BP)
set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -flto -g -fPIC")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(LINKER_OPTIONS                  -flto -Wl,--no-as-needed)

find_package(Gurobi)
find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
include_directories(${SCIP_INCLUDE_DIRS})
FILE(GLOB SOURCE_FILES src/*.cpp src/*.c src/*/*.c src/*/*.cpp src/*/*/*.c src/*/*/*.cpp)
add_executable(BP ${SOURCE_FILES})

target_link_libraries(  BP PRIVATE ${GUROBI_CXX_LIBRARY})
target_link_libraries(  BP PRIVATE ${GUROBI_LIBRARY})
target_link_libraries(  BP PRIVATE dl   
			${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(BP PRIVATE scip)