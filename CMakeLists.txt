cmake_minimum_required(VERSION 2.6)
project(ASTAR-JSP)

set (SOURCES AStar.cpp IndexPriorityQueue.cpp)
include_directories(${ASTAR-JSP_SOURCE_DIR})

add_library (libastarjsp ${SOURCES})

add_executable (astar_jsp_test TestAStar.cpp)

target_link_libraries (astar_jsp_test libastarjsp)

