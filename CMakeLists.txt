cmake_minimum_required(VERSION 2.8.11)
project(BlendshapeGeneration)

set(CUDA_USE_STATIC_CUDA_RUNTIME OFF)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") #c++11
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

message(STATUS "The CXX flags: ${CMAKE_CXX_FLAGS}")

# OpenMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories("/usr/local/include")


# Eigen
set(EIGEN3_INCLUDE_DIRS /usr/local/include/eigen3/)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message("--   Eigen include path is ${EIGEN3_INCLUDE_DIRS}")

# SuiteSparse
set(SUITESPARSE_INCLUDE_DIR "/usr/include/suitesparse")
include_directories(${SUITESPARSE_INCLUDE_DIR})

find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
if(NOT SuiteSparse_FOUND)
        #set(SuiteSparse_VERBOSE ON)
        find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
else()
        message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
        include(${USE_SuiteSparse})
endif()
MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")

# Ceres solver
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
link_libraries(${CERES_LIBRARIES})

# trimesh
set(Trimesh2_DIR ${PROJECT_SOURCE_DIR}/third_party/trimesh2)
set(Trimesh2_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/third_party/trimesh2/include/)
set(Trimesh2_LIB_DIRS ${PROJECT_SOURCE_DIR}/third_party/trimesh2/lib.Linux64/libtrimesh.a)
include_directories(${Trimesh2_INCLUDE_DIRS})

# json
include_directories("${PROJECT_SOURCE_DIR}/third_party/json/include")


###
include_directories(${PROJECT_SOURCE_DIR}/include)

set(ALL_NEED_LIBS
${CERES_LIBRARIES}
${SuiteSparse_LIBRARIES}
${Trimesh2_LIB_DIRS} 
)

add_library(${PROJECT_NAME} STATIC 
src/TriangleGradient.cpp
src/MeshTransfer.cpp 
src/BlendshapeGenerator.cpp
)

target_link_libraries(${PROJECT_NAME} ${ALL_NEED_LIBS})


add_executable(TestCases src/TestCases.cpp)
target_link_libraries(TestCases ${PROJECT_NAME})
