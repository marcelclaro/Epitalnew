cmake_minimum_required(VERSION 2.8)
set (project epital)

set(${project}_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PerturbationMethod.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverTM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PWFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverSO3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Materials.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CMaterial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Heterostructure.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Units.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Heterostructure3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverSO.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WannierFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Graphics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Alloy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lab3D.cpp
)

set(${project}_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Materials.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Units.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solution3D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolutionWannier.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WannierFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PerturbationMethod.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Graphics.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CMaterial.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Constants.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction3D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction2D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverTM.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution3D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution2D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CustomTypes.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolutionPW.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverSO3D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Grid.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Heterostructure.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Heterostructure3D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PWFunction.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WaveTimeEvolution.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solution.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SolverSO.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BrentsRootFinder.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Alloy.hpp
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_mod/")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK REQUIRED)
find_package(FFTW REQUIRED)
find_package(OpenMP)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_definitions( "-DARMA_USE_HDF5" )
include(${VTK_USE_FILE})

add_executable(${project} ${${project}_SOURCES} ${${project}_HEADERS} )
target_link_libraries (${project} Eigen3::Eigen)
target_link_libraries (${project} ${VTK_LIBRARIES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(${project} PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(${project} ${ARMADILLO_LIBRARIES})
