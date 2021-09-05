cmake_minimum_required(VERSION 3.1)

project(CG)

set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -flto -g -fPIC -fopenmp")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(LINKER_OPTIONS                  -flto -Wl,--no-as-needed)

FILE(GLOB SOURCE_FILES src/*.cpp src/*.cc src/*/*.c src/*/*.cpp src/*/*/*.c)

find_package(Gurobi)
include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
include_directories(SYSTEM "/usr/local/lib/")

add_executable(CG ${SOURCE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/prebuild/fastcluster.o)

target_include_directories(CG PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/svm"
				       "${CMAKE_CURRENT_SOURCE_DIR}/src/linear_svm"
				       "${CMAKE_CURRENT_SOURCE_DIR}/src/linear_svm/blas")
target_link_libraries(  CG PRIVATE ${GUROBI_CXX_LIBRARY})
target_link_libraries(  CG PRIVATE ${GUROBI_LIBRARY})
target_link_libraries(  CG PRIVATE dl)

find_package(Boost COMPONENTS system filesystem REQUIRED)
target_link_libraries(  CG PRIVATE dl   
			${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
