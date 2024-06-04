cmake_minimum_required(VERSION 3.27)

project(RedimesionnementPixel)

cmake_policy(SET CMP0057 NEW)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_definitions(-g)

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES})

#add_executable(Algo1 algo1.cpp CLI11.hpp DoublyLinkedList.h DoublyLinkedList.cpp)
#target_link_libraries(Algo1 ${DGTAL_LIBRARIES})
