#include <gtest/gtest.h>
#include <mange/SO2.h>

#include <eigen3/Eigen/Dense>

#include <random>
#include <vector>

// #define ASSERT_MATRIX_EQ(lhs, rhs) \
//   for (size_t i = 0; i < (lhs).rows(); i++) { \
//     for (size_t j = 0; j < (lhs).cols(); j++) { \
//       ASSERT_DOUBLE_EQ((lhs)(i,j), (rhs)(i,j)); \
//     } \
//   }

#define ASSERT_MATRIX_EQ(actual, expect) ASSERT_TRUE((actual).isApprox(expect))
#define EXPECT_MATRIX_EQ(actual, expect) EXPECT_TRUE((actual).isApprox(expect))

#define ASSERT_SO2_EQ(actual, expect) ASSERT_MATRIX_EQ((actual).C(), (expect).C())

#define ASSERT_MATRIX_ZERO(actual) ASSERT_TRUE((actual).isMuchSmallerThan())
#define EXPECT_MATRIX_ZERO(actual) EXPECT_TRUE((actual).isMuchSmallerThan())

double wrap_angle(double angle)
{
  while (angle >  M_PI) angle -= 2*M_PI;
  while (angle < -M_PI) angle += 2*M_PI;
  return angle;
}

class SO2Test : public ::testing::Test
{
protected:
  static constexpr size_t SIZE = 100;

  std::vector<double> random_phi;
  std::vector<mange::SO2> random_X;

  SO2Test()
  {
    for (size_t i = 0; i < SIZE; i++)
    {
      std::random_device rd;
      std::mt19937_64 gen(rd());
      std::uniform_real_distribution<> dist(-100.0, 100.0);

      random_phi.push_back(dist(gen));
      random_X.emplace_back(dist(gen));
    }
  }
};

TEST_F(SO2Test, DefaultValue)
{
  mange::SO2 X;
  ASSERT_MATRIX_EQ(X.C(), Eigen::Matrix2d::Identity());
}

TEST_F(SO2Test, IdentityValue)
{
  ASSERT_MATRIX_EQ(mange::SO2::Identity().C(), Eigen::Matrix2d::Identity());
}

TEST_F(SO2Test, SpecialOrthogonal)
{
  for (const auto& X : random_X)
  {
    ASSERT_MATRIX_EQ(X.C() * X.C().transpose(), Eigen::Matrix2d::Identity());
    ASSERT_MATRIX_EQ(X.C().transpose() * X.C(), Eigen::Matrix2d::Identity());
    ASSERT_DOUBLE_EQ(X.C().determinant(), 1.0);
  }
}

TEST_F(SO2Test, Closure)
{
  for (size_t i = 0; i < random_X.size(); i++)
  {
    mange::SO2 X = random_X[i] * random_X[random_X.size() - i - 1];

    // make sure it's still a special orthogonal matrix
    ASSERT_MATRIX_EQ(X.C() * X.C().transpose(), Eigen::Matrix2d::Identity());
    ASSERT_MATRIX_EQ(X.C().transpose() * X.C(), Eigen::Matrix2d::Identity());
    ASSERT_DOUBLE_EQ(X.C().determinant(), 1.0);
  }
}

TEST_F(SO2Test, Identity)
{
  mange::SO2 identity;
  for (const auto& X : random_X)
  {
    ASSERT_SO2_EQ(X * identity, X);
    ASSERT_SO2_EQ(identity * X, X);
  }
}

TEST_F(SO2Test, Inverse)
{
  for (const auto& X : random_X)
  {
    ASSERT_MATRIX_EQ((X * X.inverse()).C(), Eigen::Matrix2d::Identity());
    ASSERT_MATRIX_EQ((X.inverse() * X).C(), Eigen::Matrix2d::Identity());
  }
}

TEST_F(SO2Test, Associativity)
{
  for (size_t i = 0; i < random_X.size() - 2; i++)
  {
    ASSERT_MATRIX_EQ(((random_X[i] * random_X[i + 1]) * random_X[i + 2]).C(),
                     (random_X[i] * (random_X[i + 1] * random_X[i + 2])).C());
  }
}

TEST_F(SO2Test, HatVeeInverseMappings)
{
  for (double phi : random_phi)
  {
    ASSERT_DOUBLE_EQ(mange::SO2::vee(mange::SO2::hat(phi)), phi);
  }

  //! @todo test other direction (would require having log() -> Matrix2d function)
}

TEST_F(SO2Test, ExpLogInverseMappings)
{
  for (double phi : random_phi)
  {
    ASSERT_NEAR(mange::SO2::Exp(phi).Log(), wrap_angle(phi), 1e-12);
  }

  for (const auto& X : random_X)
  {
    ASSERT_SO2_EQ(mange::SO2::Exp(X.Log()), X);
  }
}

TEST_F(SO2Test, Adjoint)
{
  for (size_t i = 0; i < SIZE; i++)
  {
    ASSERT_DOUBLE_EQ(random_X[i].Ad() * random_phi[i],
                     mange::SO2::vee(random_X[i].C() * mange::SO2::hat(random_phi[i]) * random_X[i].inverse().C()));
  }
}