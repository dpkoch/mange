#include "test/helpers.hh"

namespace test {

bool doubleEqual(double lhs, double rhs) {
    ::testing::internal::Double left(lhs), right(rhs);
    return left.AlmostEquals(right);
}

bool doubleNear(double lhs, double rhs, double epsilon) {
    return std::abs(lhs - rhs) <= std::max(std::abs(rhs), std::abs(lhs)) * epsilon;
}

::testing::AssertionResult nearlyEqual(double lhs, double rhs, double epsilon) {
    if (doubleNear(lhs, rhs, epsilon))
        return ::testing::AssertionSuccess();
    else {
        return ::testing::AssertionFailure() << "lhs: " << lhs << std::endl << "rhs: " << rhs;
    }
}

}  // namespace test