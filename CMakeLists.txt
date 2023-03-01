cmake_minimum_required(VERSION 3.20)


set(This Heisenberg)

project(${This} VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)






add_subdirectory(HighFive)




add_executable(${This} main.cpp)
target_link_libraries(${This} HighFive)


find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${This} OpenMP::OpenMP_CXX)
endif()



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
