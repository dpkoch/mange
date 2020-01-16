#include <gtest/gtest.h>
#include <mange/SO2.h>

#include <eigen3/Eigen/Dense>

#include <cmath>
#include <random>
#include <vector>

//==============================================================================
// type definitions
//==============================================================================

// This struct defines three types associated with a group Lie group:
//   Vector: the type of the vector space isomorphic to the Lie algebra
//   Algebra: the type of the Lie algebra
//   Domain: the type of the domain group action (e.g. rotation or translation of vectors)
template <typename Group> struct TypeDefs;

template <> struct TypeDefs<mange::SO2>
{
  typedef double Vector;
  typedef Eigen::Matrix2d Algebra;
  typedef Eigen::Vector2d Domain;
};

//==============================================================================
// helper functions
//==============================================================================

bool doubleEqual(double lhs, double rhs)
{
  ::testing::internal::Double left(lhs), right(rhs);
  return left.AlmostEquals(right);
}

bool doubleNear(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon())
{
  return std::abs(lhs - rhs) <= ( std::abs(lhs) > std::abs(rhs) ? std::abs(lhs) : std::abs(rhs) ) * epsilon;
}

// valid group member
template <typename Group>
bool isValidGroupMember(const Group &X) = delete;

template <>
bool isValidGroupMember(const mange::SO2 &X)
{
  // should be orthogonal matrix with determinant 1
  return (X.C() * X.C().transpose()).isIdentity()
         && (X.C().transpose() * X.C()).isIdentity()
         && doubleEqual(X.C().determinant(), 1.0);
}

// identity element
template <typename Group>
bool isIdentityElement(const Group &X) = delete;

template <>
bool isIdentityElement(const mange::SO2 &X)
{
  return X.C().isIdentity();
}

// random vectors
template <typename Vector>
Vector randomVector()
{
  return Vector::Random();
}

template <>
double randomVector()
{
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<> dist(-100.0, 100.0);

  return dist(gen);
}

template <typename Domain>
auto constexpr randomDomain = &randomVector<Domain>;

// equality of vectors
template <typename Vector>
bool vectorsEqual(const Vector &actual, const Vector &expect)
{
  return actual.isApprox(expect);
}

bool vectorsEqual(double actual, double expect)
{
  return doubleEqual(actual, expect);
}

// compare constrained tangent space elements
//   These functions are used to deal with the many-to-one property of the
//   exponential mapping. As a result of this property Log(Exp(x)) will not
//   necessarily return x. These functions map an arbitrary vector x to an
//   equivalent vector in the range of Log().

double wrap_angle(double angle)
{
  // add/subtract multiples of 2pi all in one step, as opposed to sequentially, to reduce error
  size_t add = 0;
  size_t subtract = 0;

  while (angle - subtract*2*M_PI > M_PI)
    ++subtract;

  while (angle + add*2*M_PI < -M_PI)
    ++add;

  return angle + add*2*M_PI - subtract*2*M_PI;
}

template <typename Group>
::testing::AssertionResult tangentVectorsEqual(const typename TypeDefs<Group>::Vector &log_x,
                                               const typename TypeDefs<Group>::Vector &original_x);

template <>
::testing::AssertionResult tangentVectorsEqual<mange::SO2>(const double &log_x, const double &original_x) // const double& is innefficient here, but required to match template
{
  double constrained_x = wrap_angle(original_x);

  if (doubleNear(log_x, constrained_x, 1e-12))
    return ::testing::AssertionSuccess();
  else
    return ::testing::AssertionFailure() << "log: " << log_x << ", constrained: " << constrained_x;
}

// check that group action does not change norm of vectors in its domain
template <typename Domain>
::testing::AssertionResult normEqual(const Domain &actual, const Domain &expect)
{
  if (doubleEqual(actual.norm(), expect.norm()))
    return ::testing::AssertionSuccess();
  else
    return ::testing::AssertionFailure() << "actual: " << actual.norm() << ", expected: " << expect.norm();
}

//==============================================================================
// test fixture
//==============================================================================

template <typename Group>
class LieGroupTest : public ::testing::Test
{
protected:
  using Vector = typename TypeDefs<Group>::Vector;
  using Domain = typename TypeDefs<Group>::Domain;

  static constexpr size_t SIZE = 100;

  std::vector<Group> random_X;
  std::vector<Vector> random_x;
  std::vector<Domain> random_domain;

  LieGroupTest()
  {
    for (size_t i = 0; i < SIZE; i++)
    {
      random_X.push_back(Group::Random());
      random_x.push_back(randomVector<Vector>());
      random_domain.push_back(randomDomain<Domain>());
    }
  }
};

//==============================================================================
// test cases
//==============================================================================

using LieGroupTypes = ::testing::Types<mange::SO2>;
TYPED_TEST_CASE(LieGroupTest, LieGroupTypes); //! @note TYPED_TEST_CASE will be replace by TYPED_TEST_SUITE in googletest v1.10.x

TYPED_TEST(LieGroupTest, IdentityValue)
{
  ASSERT_TRUE(isIdentityElement(TypeParam::Identity()));
}

TYPED_TEST(LieGroupTest, DefaultValue)
{
  ASSERT_TRUE(isIdentityElement(TypeParam()));
}

TYPED_TEST(LieGroupTest, Closure)
{
  for (size_t i = 0; i < this->random_X.size(); i++)
  {
    TypeParam X = this->random_X[i] * this->random_X[this->random_X.size() - i - 1];
    ASSERT_TRUE(isValidGroupMember(X));
  }
}

TYPED_TEST(LieGroupTest, Identity)
{
  TypeParam identity;
  for (const auto &X : this->random_X)
  {
    ASSERT_TRUE((X * identity).isApprox(X));
    ASSERT_TRUE((identity * X).isApprox(X));
  }
}

TYPED_TEST(LieGroupTest, Inverse)
{
  for (const auto &X : this->random_X)
  {
    ASSERT_TRUE((X * X.inverse()).isIdentity());
    ASSERT_TRUE((X.inverse() * X).isIdentity());
  }
}

TYPED_TEST(LieGroupTest, Associativity)
{
  for (size_t i = 0; i < this->random_X.size() - 2; i++)
  {
    const auto& A = this->random_X[i];
    const auto& B = this->random_X[i+1];
    const auto& C = this->random_X[i+2];

    ASSERT_TRUE(((A * B) * C).isApprox(A * (B * C)));
  }
}

TYPED_TEST(LieGroupTest, HatVeeInverseMappings)
{
  constexpr auto hat = &TypeParam::hat;
  constexpr auto vee = &TypeParam::vee;

  for (const auto &x : this->random_x)
  {
    ASSERT_TRUE(vectorsEqual(vee(hat(x)), x));
  }

  //! @todo test hat(vee(x)) direction (would require having log() -> Algebra function)
}

TYPED_TEST(LieGroupTest, ExpLogInverseMappings)
{
  constexpr auto Exp = &TypeParam::Exp;

  for (const auto &x : this->random_x)
  {
    ASSERT_TRUE(tangentVectorsEqual<TypeParam>(Exp(x).Log(), x));
  }

  for (const auto &X : this->random_X)
  {
    ASSERT_TRUE(Exp(X.Log()).isApprox(X));
  }
}

TYPED_TEST(LieGroupTest, Adjoint)
{
  constexpr auto hat = &TypeParam::hat;
  constexpr auto vee = &TypeParam::vee;

  for (size_t i = 0; i < TestFixture::SIZE; i++)
  {
    const auto &X = this->random_X[i];
    const auto &x = this->random_x[i];

    ASSERT_TRUE(vectorsEqual(X.Ad() * x, vee(X.matrix() * hat(x) * X.inverse().matrix())));
  }
}

TYPED_TEST(LieGroupTest, Action)
{
  for (size_t i = 0; i < TestFixture::SIZE; i++)
  {
    const auto &X = this->random_X[i];
    const auto &v = this->random_domain[i];

    ASSERT_TRUE(normEqual(X*v, v));
  }
}