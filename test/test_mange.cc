#include <cmath>
#include <random>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <gtest/gtest.h>

#include "mange/mange.hh"

//==============================================================================
// helper functions
//==============================================================================

bool doubleEqual(double lhs, double rhs) {
    ::testing::internal::Double left(lhs), right(rhs);
    return left.AlmostEquals(right);
}

bool doubleNear(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon()) {
    return std::abs(lhs - rhs) <= std::max(std::abs(rhs), std::abs(lhs)) * epsilon;
}

template <typename DerivedLhs, typename DerivedRhs>
::testing::AssertionResult nearlyEqual(const Eigen::MatrixBase<DerivedLhs> &lhs,
                                       const Eigen::MatrixBase<DerivedRhs> &rhs,
                                       double epsilon = 1e-12) {
    if (lhs.isApprox(rhs))
        return ::testing::AssertionSuccess();
    else {
        return ::testing::AssertionFailure() << "lhs:" << std::endl
                                             << lhs << std::endl
                                             << "rhs:" << std::endl
                                             << rhs;
    }
}

::testing::AssertionResult nearlyEqual(double lhs, double rhs, double epsilon = 1e-12) {
    if (doubleNear(lhs, rhs, epsilon))
        return ::testing::AssertionSuccess();
    else {
        return ::testing::AssertionFailure() << "lhs: " << lhs << std::endl << "rhs: " << rhs;
    }
}

// valid group member
template <typename Group>
::testing::AssertionResult isValidSOnElement(const Group &X) {
    const auto &matrix = X.matrix();

    if (!(matrix * matrix.transpose()).isIdentity()) {
        return ::testing::AssertionFailure() << "C * C^T is not identity: " << std::endl
                                             << matrix * matrix.transpose();
    }
    if (!(matrix.transpose() * matrix).isIdentity()) {
        return ::testing::AssertionFailure() << "C^T * C is not identity: " << std::endl
                                             << matrix.transpose() * matrix;
    }
    if (!doubleEqual(matrix.determinant(), 1.0)) {
        return ::testing::AssertionFailure() << "Determinant is not 1: " << matrix.determinant();
    }

    return ::testing::AssertionSuccess();
}

template <typename Group>
::testing::AssertionResult isValidSEnElement(const Group &X) {
    ::testing::AssertionResult rotation_result = isValidSOnElement(X.C());
    if (!rotation_result) return rotation_result;

    const auto matrix = X.matrix();
    if (!matrix.template topLeftCorner<Group::DIM, Group::DIM>().isApprox(X.C().matrix())) {
        return ::testing::AssertionFailure()
               << "Top left corner not equal to rotation matrix" << std::endl
               << "Top left corner:" << std::endl
               << matrix.template topLeftCorner<Group::DIM, Group::DIM>() << std::endl
               << "Rotation matrix:" << std::endl
               << X.C().matrix();
    }
    if (!matrix.template bottomLeftCorner<1, Group::DIM>().isZero()) {
        return ::testing::AssertionFailure() << "Bottom left corner is not zero:" << std::endl
                                             << matrix.template bottomLeftCorner<1, Group::DIM>();
    }
    if (!(matrix.template bottomRightCorner<1, 1>()(0, 0) == 1.0)) {
        return ::testing::AssertionFailure() << "Bottom right corner is not 1: "
                                             << matrix.template bottomRightCorner<1, 1>()(0, 0);
    }

    return ::testing::AssertionSuccess();
}

template <typename Group>
::testing::AssertionResult isValidGroupMember(const Group &X) = delete;

template <>
::testing::AssertionResult isValidGroupMember(const mange::SO2 &X) {
    return isValidSOnElement(X);
}

template <>
::testing::AssertionResult isValidGroupMember(const mange::SO3 &X) {
    return isValidSOnElement(X);
}

template <>
::testing::AssertionResult isValidGroupMember(const mange::SE2 &X) {
    return isValidSEnElement(X);
}

template <>
::testing::AssertionResult isValidGroupMember(const mange::SE3 &X) {
    return isValidSEnElement(X);
}

// identity element
template <typename Group>
bool isIdentityElement(const Group &X) = delete;

template <>
bool isIdentityElement(const mange::SO2 &X) {
    return X.matrix().isIdentity();
}

template <>
bool isIdentityElement(const mange::SE2 &X) {
    return isIdentityElement(X.C()) && X.r().isZero();
}

template <>
bool isIdentityElement(const mange::SO3 &X) {
    return X.matrix().isIdentity();
}

template <>
bool isIdentityElement(const mange::SE3 &X) {
    return isIdentityElement(X.C()) && X.r().isZero();
}

// zero vectors
template <typename Vector>
Vector zero_vector() {
    return Vector::Zero();
}

template <>
mange::SO2::VectorType zero_vector<mange::SO2::VectorType>() {
    return 0.0;
}

// random vectors
template <typename Vector>
Vector randomVector() {
    return Vector::Random();
}

template <>
double randomVector() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dist(-100.0, 100.0);

    return dist(gen);
}

template <typename Domain>
auto constexpr randomDomain = &randomVector<Domain>;

// equality of vectors
template <typename Vector>
bool vectorsEqual(const Vector &actual, const Vector &expect) {
    return actual.isApprox(expect);
}

bool vectorsEqual(double actual, double expect) {
    return doubleEqual(actual, expect);
}

// compare constrained tangent space elements
//   These functions are used to deal with the many-to-one property of the
//   exponential mapping. As a result of this property Log(Exp(x)) will not
//   necessarily return x. These functions map an arbitrary vector x to an
//   equivalent vector in the range of Log().

double wrap_angle(double angle) {
    // add/subtract multiples of 2pi all in one step, as opposed to sequentially, to reduce error
    int add = 0;
    int subtract = 0;

    while (angle - subtract * 2 * M_PI > M_PI) ++subtract;

    while (angle + add * 2 * M_PI < -M_PI) ++add;

    return angle + (add - subtract) * 2 * M_PI;
}

mange::SO3::VectorType wrap_angle(const mange::SO3::VectorType &phi) {
    return wrap_angle(phi.norm()) * phi.normalized();
}

template <typename Group>
::testing::AssertionResult tangentVectorsEqual(const typename Group::VectorType &log_x,
                                               const typename Group::VectorType &original_x);

template <>
::testing::AssertionResult tangentVectorsEqual<mange::SO2>(
    const double &log_x,
    const double &original_x)  // const double& is innefficient here, but required to match template
{
    double constrained_x = wrap_angle(original_x);

    if (doubleNear(log_x, constrained_x, 1e-10))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure() << "log: " << std::endl
                                             << log_x << ", constrained: " << std::endl
                                             << constrained_x;
}

template <>
::testing::AssertionResult tangentVectorsEqual<mange::SE2>(
    const mange::SE2::VectorType &log_x, const mange::SE2::VectorType &original_x) {
    mange::SE2::VectorType constrained_x;
    constrained_x << original_x(0), original_x(1), wrap_angle(original_x(2));

    if (original_x.isApprox(log_x))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure() << "log: " << std::endl
                                             << log_x << ", constrained: " << std::endl
                                             << constrained_x;
}

template <>
::testing::AssertionResult tangentVectorsEqual<mange::SO3>(
    const mange::SO3::VectorType &log_x, const mange::SO3::VectorType &original_x) {

    mange::SO3::VectorType constrained_x = wrap_angle(original_x);

    if (nearlyEqual(log_x, constrained_x, 1e-10))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure()
               << "log: " << log_x << ", constrained: " << constrained_x;
}

template <>
::testing::AssertionResult tangentVectorsEqual<mange::SE3>(
    const mange::SE3::VectorType &log_x, const mange::SE3::VectorType &original_x) {
    mange::SE3::VectorType constrained_x = original_x;
    mange::SE3::angular_part(constrained_x) = wrap_angle(mange::SE3::angular_part(original_x));

    if (original_x.isApprox(log_x))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure()
               << "log: " << log_x << ", constrained: " << constrained_x;
}

// check that group action satisfies any constraints

template <typename Domain>
::testing::AssertionResult normEqual(const Domain &after, const Domain &before) {
    if (doubleEqual(after.norm(), before.norm()))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure()
               << "after: " << after.norm() << ", before: " << before.norm();
}

template <typename Group>
::testing::AssertionResult actionValid(const typename Group::DomainType &after,
                                       const typename Group::DomainType &before) {
    return ::testing::AssertionSuccess();  // no checks for SE(n)
}

template <>
::testing::AssertionResult actionValid<mange::SO2>(const mange::SO2::DomainType &after,
                                                   const mange::SO2::DomainType &before) {
    return normEqual(after, before);
}

template <>
::testing::AssertionResult actionValid<mange::SO3>(const mange::SO3::DomainType &after,
                                                   const mange::SO3::DomainType &before) {
    return normEqual(after, before);
}

// check that (left or right) Jacobian and its inverse are actually inverses
template <typename T>
::testing::AssertionResult jacobianAndInverse(const T &J, const T &Jinverse) {
    if ((J * Jinverse).isIdentity(1e-10) && (Jinverse * J).isIdentity(1e-10))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure()
               << "J and Jinverse are not inverses of each other!" << std::endl
               << "J: " << std::endl
               << J << std::endl
               << "Jinverse: " << std::endl
               << Jinverse << std::endl
               << "J*Jinverse: " << std::endl
               << (J * Jinverse) << std::endl
               << "Jinverse*J: " << std::endl
               << (Jinverse * J);
}

::testing::AssertionResult jacobianAndInverse(double J, double Jinverse) {
    if (doubleEqual(J * Jinverse, 1.0))
        return ::testing::AssertionSuccess();
    else
        return ::testing::AssertionFailure()
               << "J and Jinverse are not inverses of each other! J: " << J
               << ", Jinverse: " << Jinverse;
}

//==============================================================================
// test fixture
//==============================================================================

enum class DataSet { ZERO, RANDOM };

template <typename G, DataSet D>
struct TestTraits {
    using Group = G;
    static constexpr DataSet DATASET = D;
};

template <typename Traits>
class LieGroupTest : public ::testing::Test {
   protected:
    using Group = typename Traits::Group;
    static constexpr DataSet DATASET = Traits::DATASET;

    using Vector = typename Group::VectorType;
    using Domain = typename Group::DomainType;

    static constexpr size_t SIZE = 100;

    std::vector<Group> random_X;
    std::vector<Vector> random_x;
    std::vector<Domain> random_domain;

    LieGroupTest() {
        switch (DATASET) {
            case DataSet::ZERO:
                random_X.push_back(Group::Identity());
                random_x.push_back(zero_vector<Vector>());
                random_domain.push_back(Domain::Zero());

                random_X.push_back(Group::Identity());
                random_x.push_back(zero_vector<Vector>());
                random_domain.push_back(Domain::Ones());
                break;
            case DataSet::RANDOM:
                random_X.reserve(SIZE);
                random_x.reserve(SIZE);
                random_domain.reserve(SIZE);
                for (size_t i = 0; i < SIZE; i++) {
                    random_X.push_back(Group::Random());
                    random_x.push_back(randomVector<Vector>());
                    random_domain.push_back(randomDomain<Domain>());
                }
                break;
        }
    }
};

//==============================================================================
// test cases
//==============================================================================

using TestTraitsList = ::testing::Types<TestTraits<mange::SO2, DataSet::ZERO>,
                                        TestTraits<mange::SO2, DataSet::RANDOM>,
                                        TestTraits<mange::SE2, DataSet::ZERO>,
                                        TestTraits<mange::SE2, DataSet::RANDOM>,
                                        TestTraits<mange::SO3, DataSet::ZERO>,
                                        TestTraits<mange::SO3, DataSet::RANDOM>,
                                        TestTraits<mange::SE3, DataSet::ZERO>,
                                        TestTraits<mange::SE3, DataSet::RANDOM>>;
TYPED_TEST_SUITE(LieGroupTest, TestTraitsList);

TYPED_TEST(LieGroupTest, IdentityValue) {
    ASSERT_TRUE(isIdentityElement(TypeParam::Group::Identity()));
}

TYPED_TEST(LieGroupTest, DefaultValue) {
    ASSERT_TRUE(isIdentityElement(typename TypeParam::Group()));
}

TYPED_TEST(LieGroupTest, Exp) {
    for (const auto &x : this->random_x) {
        typename TypeParam::Group X = TypeParam::Group::Exp(x);
        ASSERT_TRUE(isValidGroupMember(X));
    }
}

TYPED_TEST(LieGroupTest, Closure) {
    for (size_t i = 0; i < this->random_X.size(); i++) {
        typename TypeParam::Group X =
            this->random_X[i] * this->random_X[this->random_X.size() - i - 1];
        ASSERT_TRUE(isValidGroupMember(X));
    }
}

TYPED_TEST(LieGroupTest, Identity) {
    typename TypeParam::Group identity;
    for (const auto &X : this->random_X) {
        ASSERT_TRUE((X * identity).isApprox(X));
        ASSERT_TRUE((identity * X).isApprox(X));
    }
}

TYPED_TEST(LieGroupTest, Inverse) {
    for (const auto &X : this->random_X) {
        ASSERT_TRUE((X * X.inverse()).isIdentity());
        ASSERT_TRUE((X.inverse() * X).isIdentity());
    }
}

TYPED_TEST(LieGroupTest, Associativity) {
    for (size_t i = 0; i < this->random_X.size() - 2; i++) {
        const auto &A = this->random_X[i];
        const auto &B = this->random_X[i + 1];
        const auto &C = this->random_X[i + 2];

        ASSERT_TRUE(((A * B) * C).isApprox(A * (B * C)));
    }
}

TYPED_TEST(LieGroupTest, HatVeeInverseMappings) {
    // using TypeParam::Group::typename hat;
    // using TypeParam::Group::typename vee;

    for (const auto &x : this->random_x) {
        ASSERT_TRUE(vectorsEqual(TypeParam::Group::vee(TypeParam::Group::hat(x)), x));
    }

    //! @todo test hat(vee(x)) direction (would require having log() -> Algebra function)
}

TYPED_TEST(LieGroupTest, ExpLogInverseMappings) {
    constexpr auto Exp = &TypeParam::Group::Exp;

    for (const auto &x : this->random_x) {
        ASSERT_TRUE(tangentVectorsEqual<typename TypeParam::Group>(Exp(x).Log(), x));
    }

    for (const auto &X : this->random_X) {
        ASSERT_TRUE(Exp(X.Log()).isApprox(X));
    }
}

TYPED_TEST(LieGroupTest, Adjoint) {
    constexpr auto hat = &TypeParam::Group::hat;
    constexpr auto vee = &TypeParam::Group::vee;

    for (size_t i = 0; i < TestFixture::SIZE; i++) {
        const auto &X = this->random_X[i];
        const auto &x = this->random_x[i];

        typename TypeParam::Group::VectorType Ad_x =
            X.Ad() * x;  // assignment required to evaluate Eigen product expression down to matrix
                         // type for template matching
        ASSERT_TRUE(vectorsEqual(Ad_x, vee(X.matrix() * hat(x) * X.inverse().matrix())));
    }
}

TYPED_TEST(LieGroupTest, Action) {
    for (size_t i = 0; i < TestFixture::SIZE; i++) {
        const auto &X = this->random_X[i];
        const auto &v = this->random_domain[i];

        ASSERT_TRUE(actionValid<typename TypeParam::Group>(X * v, v));
    }
}

TYPED_TEST(LieGroupTest, IdentityAction) {
    const auto X = TypeParam::Group::Identity();
    for (size_t i = 0; i < TestFixture::SIZE; i++) {
        const auto &v = this->random_domain[i];

        ASSERT_TRUE(nearlyEqual(X * v, v));
    }
}

TYPED_TEST(LieGroupTest, JacobiansAndInverses) {
    for (const auto &x : this->random_x) {
        ASSERT_TRUE(jacobianAndInverse(TypeParam::Group::Jl(x), TypeParam::Group::JlInverse(x)));
        ASSERT_TRUE(jacobianAndInverse(TypeParam::Group::Jr(x), TypeParam::Group::JrInverse(x)));
    }
}

TYPED_TEST(LieGroupTest, AdjointAndJacobians) {
    for (const auto &x : this->random_x) {
        const typename TypeParam::Group X = TypeParam::Group::Exp(x);
        ASSERT_TRUE(nearlyEqual(X.Ad(), TypeParam::Group::Jl(x) * TypeParam::Group::JrInverse(x)));
    }
}

TYPED_TEST(LieGroupTest, InverseOfExpIsExpOfNegative) {
    for (const auto &x : this->random_x) {
        ASSERT_TRUE(TypeParam::Group::Exp(x).inverse().isApprox(TypeParam::Group::Exp(-x)));
    }
}