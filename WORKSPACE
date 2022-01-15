load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

GOOGLETEST_TAG = "release-1.11.0"

GOOGLETEST_SHA = "353571c2440176ded91c2de6d6cd88ddd41401d14692ec1f99e35d013feda55a"

http_archive(
    name = "googletest",
    sha256 = GOOGLETEST_SHA,
    strip_prefix = "googletest-{}".format(GOOGLETEST_TAG),
    urls = ["https://github.com/google/googletest/archive/refs/tags/{}.zip".format(GOOGLETEST_TAG)],
)

EIGEN_TAG = "3.4.0"

EIGEN_SHA = "1ccaabbfe870f60af3d6a519c53e09f3dcf630207321dffa553564a8e75c4fc8"

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
