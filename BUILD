load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "mange",
    srcs = [
        "src/mange/SE2.cc",
        "src/mange/SE3.cc",
        "src/mange/SO2.cc",
        "src/mange/SO3.cc",
    ],
    hdrs = [
        "include/mange/SE2.hh",
        "include/mange/SE3.hh",
        "include/mange/SO2.hh",
        "include/mange/SO3.hh",
        "include/mange/mange.hh",
    ],
    strip_include_prefix = "include",
    visibility = ["//visibility:public"],
    deps = ["@eigen3"],
)

cc_test(
    name = "test_mange",
    size = "small",
    srcs = ["test/test_mange.cc"],
    deps = [
        ":mange",
        "@googletest//:gtest_main",
    ],
)
