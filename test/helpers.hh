#include <Eigen/Dense>
#include <gtest/gtest.h>

namespace test {

bool doubleEqual(double lhs, double rhs);

bool doubleNear(double lhs, double rhs, double epsilon = std::numeric_limits<double>::epsilon());

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

::testing::AssertionResult nearlyEqual(double lhs, double rhs, double epsilon = 1e-12);

}  // namespace test