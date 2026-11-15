#include <cmath>

#include <Eigen/Dense>

#include "mange/SE3.hh"
#include "mange/SO3.hh"
#include "test/helpers.hh"

namespace {

using mange::SE3;
using mange::SO3;

const auto X = Eigen::Vector3d::UnitX();
const auto Y = Eigen::Vector3d::UnitY();
const auto Z = Eigen::Vector3d::UnitZ();

static constexpr double SQRT2_2 = M_SQRT2 / 2.0;

void test_so3_exp(const SO3::VectorType& phi,
                  const Eigen::Vector3d& initial_vector,
                  const Eigen::Vector3d& expected) {
    EXPECT_TRUE(test::nearlyEqual(SO3::Exp(phi) * initial_vector, expected));
}

void test_se3_exp_rotation(const SO3::VectorType& phi,
                           const Eigen::Vector3d& initial_vector,
                           const Eigen::Vector3d& expected) {
    SE3::VectorType xi;
    SE3::angular_part(xi) = phi;
    SE3::linear_part(xi).setZero();

    EXPECT_TRUE(test::nearlyEqual(SE3::Exp(xi) * initial_vector, expected));
}

}  // namespace

TEST(TransformsTest, SO3Rotation) {
    // Rotate around x-axis by 90deg
    test_so3_exp(M_PI_2 * X, X, X);
    test_so3_exp(M_PI_2 * X, Y, Z);
    test_so3_exp(M_PI_2 * X, Z, -Y);

    // Rotate around x-axis by 45deg
    test_so3_exp(M_PI_4 * X, X, X);
    test_so3_exp(M_PI_4 * X, Y, SQRT2_2 * (Y + Z));
    test_so3_exp(M_PI_4 * X, Z, SQRT2_2 * (Z - Y));

    // Rotate around y-axis by 90deg
    test_so3_exp(M_PI_2 * Y, X, -Z);
    test_so3_exp(M_PI_2 * Y, Y, Y);
    test_so3_exp(M_PI_2 * Y, Z, X);

    // Rotate around y-axis by 45deg
    test_so3_exp(M_PI_4 * Y, X, SQRT2_2 * (X - Z));
    test_so3_exp(M_PI_4 * Y, Y, Y);
    test_so3_exp(M_PI_4 * Y, Z, SQRT2_2 * (Z + X));

    // Rotate around z-axis by 90deg
    test_so3_exp(M_PI_2 * Z, X, Y);
    test_so3_exp(M_PI_2 * Z, Y, -X);
    test_so3_exp(M_PI_2 * Z, Z, Z);

    // Rotate around z-axis by 45deg
    test_so3_exp(M_PI_4 * Z, X, SQRT2_2 * (X + Y));
    test_so3_exp(M_PI_4 * Z, Y, SQRT2_2 * (Y - X));
    test_so3_exp(M_PI_4 * Z, Z, Z);
}

TEST(TransformsTest, SE3Rotation) {
    // Rotate around x-axis by 90deg
    test_se3_exp_rotation(M_PI_2 * X, X, X);
    test_se3_exp_rotation(M_PI_2 * X, Y, Z);
    test_se3_exp_rotation(M_PI_2 * X, Z, -Y);

    // Rotate around x-axis by 45deg
    test_se3_exp_rotation(M_PI_4 * X, X, X);
    test_se3_exp_rotation(M_PI_4 * X, Y, SQRT2_2 * (Y + Z));
    test_se3_exp_rotation(M_PI_4 * X, Z, SQRT2_2 * (Z - Y));

    // Rotate around y-axis by 90deg
    test_se3_exp_rotation(M_PI_2 * Y, X, -Z);
    test_se3_exp_rotation(M_PI_2 * Y, Y, Y);
    test_se3_exp_rotation(M_PI_2 * Y, Z, X);

    // Rotate around y-axis by 45deg
    test_se3_exp_rotation(M_PI_4 * Y, X, SQRT2_2 * (X - Z));
    test_se3_exp_rotation(M_PI_4 * Y, Y, Y);
    test_se3_exp_rotation(M_PI_4 * Y, Z, SQRT2_2 * (Z + X));

    // Rotate around z-axis by 90deg
    test_se3_exp_rotation(M_PI_2 * Z, X, Y);
    test_se3_exp_rotation(M_PI_2 * Z, Y, -X);
    test_se3_exp_rotation(M_PI_2 * Z, Z, Z);

    // Rotate around z-axis by 45deg
    test_se3_exp_rotation(M_PI_4 * Z, X, SQRT2_2 * (X + Y));
    test_se3_exp_rotation(M_PI_4 * Z, Y, SQRT2_2 * (Y - X));
    test_se3_exp_rotation(M_PI_4 * Z, Z, Z);
}