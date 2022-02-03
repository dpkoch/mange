load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

GOOGLETEST_TAG = "release-1.11.0"

GOOGLETEST_SHA = "353571c2440176ded91c2de6d6cd88ddd41401d14692ec1f99e35d013feda55a"

http_archive(
    name = "googletest",
    sha256 = GOOGLETEST_SHA,
    strip_prefix = "googletest-{}".format(GOOGLETEST_TAG),
    urls = ["https://github.com/google/googletest/archive/refs/tags/{}.zip".format(GOOGLETEST_TAG)],
)

# Eigen version 3.4.0 breaks pybind11
EIGEN_TAG = "3.3.9"

EIGEN_SHA = "83709a8def0d60dc4d17a749989893ea5e5aacf13f9184ae0509313f400f6f45"

http_archive(
    name = "eigen3",
    build_file_content = """
cc_library(
    name = 'eigen3',
    includes = ['.'],
    hdrs = glob(['Eigen/**']),
    visibility = ['//visibility:public'],
)
""",
    sha256 = EIGEN_SHA,
    strip_prefix = "eigen-{}".format(EIGEN_TAG),
    urls = ["https://gitlab.com/libeigen/eigen/-/archive/{0}/eigen-{0}.zip".format(EIGEN_TAG)],
)

PYBIND11_BAZEL_COMMIT = "72cbbf1fbc830e487e3012862b7b720001b70672"

PYBIND11_BAZEL_SHA = "fec6281e4109115c5157ca720b8fe20c8f655f773172290b03f57353c11869c2"

http_archive(
    name = "pybind11_bazel",
    sha256 = PYBIND11_BAZEL_SHA,
    strip_prefix = "pybind11_bazel-{}".format(PYBIND11_BAZEL_COMMIT),
    urls = ["https://github.com/pybind/pybind11_bazel/archive/{}.zip".format(PYBIND11_BAZEL_COMMIT)],
)

PYBIND11_VERSION = "2.9.1"

PYBIND11_SHA = "c6160321dc98e6e1184cc791fbeadd2907bb4a0ce0e447f2ea4ff8ab56550913"

http_archive(
    name = "pybind11",
    build_file = "@pybind11_bazel//:pybind11.BUILD",
    sha256 = PYBIND11_SHA,
    strip_prefix = "pybind11-{}".format(PYBIND11_VERSION),
    urls = ["https://github.com/pybind/pybind11/archive/v{}.tar.gz".format(PYBIND11_VERSION)],
)

load("@pybind11_bazel//:python_configure.bzl", "python_configure")

python_configure(
    name = "local_config_python",
    python_version = "3",
)

http_archive(
    name = "rules_python",
    sha256 = "a30abdfc7126d497a7698c29c46ea9901c6392d6ed315171a6df5ce433aa4502",
    strip_prefix = "rules_python-0.6.0",
    url = "https://github.com/bazelbuild/rules_python/archive/0.6.0.tar.gz",
)

load("@rules_python//python:pip.bzl", "pip_install")

pip_install(
    requirements = "//:requirements.txt",
)
