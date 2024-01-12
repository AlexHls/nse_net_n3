cmake_minimum_required(VERSION 3.28)
project(nse_net_n3)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SOURCES "src/*.cpp")
include_directories(src)

option(NETWORK_SCREENING "Enable network screening" ON)
if (NETWORK_SCREENING)
    add_definitions(-DNETWORK_SCREENING)
endif()

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_executable(nse_net_n3 ${SOURCES})
target_link_libraries(nse_net_n3 ${GSL_LIBRARIES})
