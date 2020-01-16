#include <gtest/gtest.h>
#include <mange/SO2.h>

#include <eigen3/Eigen/Dense>

#include <random>
#include <vector>

::testing::AssertionResult SO2Equal(const mange::SO2 &actual, const mange::SO2 &expect)
{
  if (actual.C().isApprox(expect.C()))
    return ::testing::AssertionSuccess();
  else
    return ::testing::AssertionFailure() << "SO2 objects are not equal";

}

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
  ASSERT_TRUE(X.C().isIdentity());
}

TEST_F(SO2Test, IdentityValue)
{
  ASSERT_TRUE(mange::SO2::Identity().C().isIdentity());
}

TEST_F(SO2Test, SpecialOrthogonal)
{
  for (const auto& X : random_X)
  {
    ASSERT_TRUE((X.C() * X.C().transpose()).isIdentity());
    ASSERT_TRUE((X.C().transpose() * X.C()).isIdentity());
    ASSERT_DOUBLE_EQ(X.C().determinant(), 1.0);
  }
}

TEST_F(SO2Test, Closure)
{
  for (size_t i = 0; i < random_X.size(); i++)
  {
    mange::SO2 X = random_X[i] * random_X[random_X.size() - i - 1];

    // make sure it's still a special orthogonal matrix
    ASSERT_TRUE((X.C() * X.C().transpose()).isIdentity());
    ASSERT_TRUE((X.C().transpose() * X.C()).isIdentity());
    ASSERT_DOUBLE_EQ(X.C().determinant(), 1.0);
  }
}

TEST_F(SO2Test, Identity)
{
  mange::SO2 identity;
  for (const auto& X : random_X)
  {
    ASSERT_TRUE(SO2Equal(X * identity, X));
    ASSERT_TRUE(SO2Equal(identity * X, X));
  }
}

TEST_F(SO2Test, Inverse)
{
  for (const auto& X : random_X)
  {
    ASSERT_TRUE(SO2Equal(X * X.inverse(), mange::SO2::Identity()));
    ASSERT_TRUE(SO2Equal(X.inverse() * X, mange::SO2::Identity()));
  }
}

TEST_F(SO2Test, Associativity)
{
  for (size_t i = 0; i < random_X.size() - 2; i++)
  {
    ASSERT_TRUE(SO2Equal((random_X[i] *  random_X[i + 1]) * random_X[i + 2],
                          random_X[i] * (random_X[i + 1]  * random_X[i + 2])));
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
    ASSERT_TRUE(SO2Equal(mange::SO2::Exp(X.Log()), X));
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