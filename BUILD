load("@rules_cc//cc:defs.bzl", "cc_library")
load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
load("@rules_python//python:defs.bzl", "py_test")

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

pybind_extension(
    name = "mange_python",
    srcs = ["src/mange_python/mange_python.cc"],
    deps = [
        ":mange",
    ],
)

py_test(
    name = "test_mange_python",
    srcs = ["test/test_mange_python.py"],
    data = [":mange_python.so"],
    deps = [
        "@pip//pypi__numpy",
        "@pip//pypi__pytest",
    ],
)
