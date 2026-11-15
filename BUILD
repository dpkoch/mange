load("@rules_cc//cc:defs.bzl", "cc_library")
load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
load("@rules_python//python:pip.bzl", "compile_pip_requirements")

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
    deps = ["@eigen"],
)

pybind_extension(
    name = "mange_python",
    srcs = ["src/mange_python/mange_python.cc"],
    visibility = ["//visibility:public"],
    deps = [
        ":mange",
    ],
)

compile_pip_requirements(
    name = "requirements",
    src = "requirements.txt",
    requirements_txt = "requirements_lock.txt",
)
