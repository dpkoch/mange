cmake_minimum_required(VERSION 3.10.2)
project(mange)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${Eigen3_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
    src/${PROJECT_NAME}/SO2.cpp
    src/${PROJECT_NAME}/SE2.cpp
    src/${PROJECT_NAME}/SO3.cpp
)

include(GoogleTest)
find_package(GTest REQUIRED)
add_executable(${PROJECT_NAME}_test
    test/lie_group.cpp
)
target_link_libraries(${PROJECT_NAME}_test
    ${PROJECT_NAME}
    GTest::Main
)
gtest_discover_tests(${PROJECT_NAME}_test)