load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "mange",
    srcs = [
        "src/mange/SE2.cpp",
        "src/mange/SO2.cpp",
        "src/mange/SO3.cpp",
    ],
    hdrs = [
        "include/mange/SE2.h",
        "include/mange/SO2.h",
        "include/mange/SO3.h",
        "include/mange/mange.h",
    ],
    strip_include_prefix = "include",
    visibility = ["//visibility:public"],
    deps = ["@eigen3"],
)

cc_test(
    name = "test_mange",
    size = "small",
    srcs = ["test/lie_group.cpp"],
    deps = [
        ":mange",
        "@googletest//:gtest_main",
    ],
)
