cmake_minimum_required(VERSION 3.10.2)
project(mange)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${Eigen3_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
    src/${PROJECT_NAME}/SO2.cc
    src/${PROJECT_NAME}/SE2.cc
    src/${PROJECT_NAME}/SO3.cc
    src/${PROJECT_NAME}/SE3.cc
)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

include(GoogleTest)
find_package(GTest 1.10 REQUIRED)
add_executable(${PROJECT_NAME}_test
    test/test_mange.cc
)
target_link_libraries(${PROJECT_NAME}_test
    ${PROJECT_NAME}
    GTest::Main
)
gtest_discover_tests(${PROJECT_NAME}_test)
set_property(TARGET ${PROJECT_NAME}_test PROPERTY CXX_STANDARD 17)
