load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")
load("@rules_python//python:defs.bzl", "py_test")

cc_library(
    name = "helpers",
    srcs = ["helpers.cc"],
    hdrs = ["helpers.hh"],
    deps = [
        "@eigen",
        "@googletest//:gtest",
    ],
)

cc_test(
    name = "test_mange",
    size = "small",
    srcs = ["test_mange.cc"],
    deps = [
        ":helpers",
        "//:mange",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "test_transforms",
    size = "small",
    srcs = ["test_transforms.cc"],
    deps = [
        "helpers",
        "//:mange",
        "@googletest//:gtest_main",
    ],
)

py_test(
    name = "test_mange_python",
    size = "small",
    srcs = ["test_mange_python.py"],
    data = ["//:mange_python.so"],
    deps = [
        "@pypi//numpy",
        "@pypi//pytest",
    ],
)
