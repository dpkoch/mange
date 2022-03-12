#include <eigen3/Eigen/Dense>
#include <gtest/gtest.h>

#include "mange/algorithm/wahba.hh"
#include "mange/mange.hh"

namespace mange {

namespace {

template <typename DerivedLhs, typename DerivedRhs>
::testing::AssertionResult nearlyEqual(const Eigen::MatrixBase<DerivedLhs> &lhs,
                                       const Eigen::MatrixBase<DerivedRhs> &rhs,
                                       double epsilon = 1e-12) {
    if (lhs.isApprox(rhs, epsilon))
        return ::testing::AssertionSuccess();
    else {
        return ::testing::AssertionFailure() << "lhs:" << std::endl
                                             << lhs << std::endl
                                             << "rhs:" << std::endl
                                             << rhs;
    }
}

}  // namespace

TEST(Wahba, RecoverOriginalRotationFromSingleSample) {
    constexpr int NUM_TESTS = 100;
    constexpr double TOLERANCE = 1e-2;

    for (int i = 0; i < NUM_TESTS; i++) {
        const SO3 R = SO3::Random();
        const SO3::DomainType v = SO3::DomainType::Random();
        const SO3::DomainType w = R * v;

        const SO3 computed_R = algorithm::wahba(w, v);
        ASSERT_TRUE(nearlyEqual(computed_R.matrix(), R.matrix(), TOLERANCE));
    }
}

}  // namespace mange