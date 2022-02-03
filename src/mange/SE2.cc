#include "mange/SE2.hh"

namespace mange {

SE2::SE2() : r_(DomainType::Zero()) {}

SE2::SE2(const SO2 &rotation) : C_(rotation) {}

SE2::SE2(const DomainType &translation) : r_(translation) {}

SE2::SE2(const SO2 &rotation, const DomainType &translation) : C_(rotation), r_(translation) {}

SE2 SE2::Identity() {
    return SE2();
}

SE2 SE2::Random() {
    return SE2(SO2::Random(), DomainType::Random());
}

SE2 SE2::Exp(const VectorType &xi) {
    SE2 result;

    const SO2::VectorType phi = angular_part(xi);
    result.C_ = SO2::Exp(phi);

    double a, b;
    if (std::abs(phi) > EPSILON) {
        a = std::sin(phi) / phi;
        b = (1.0 - std::cos(phi)) / phi;
    } else  // use Taylor-series expansion
    {
        a = 1.0 - std::pow(phi, 2) / 6.0 + std::pow(phi, 4) / 120.0;
        b = phi / 2.0 - std::pow(phi, 3) / 24.0;
    }
    const Eigen::Matrix2d V = (a * Eigen::Matrix2d::Identity() + b * SO2::hat(1.0));
    result.r_ = V * linear_part(xi);

    return result;
}

SE2::VectorType SE2::Log() const {
    const SO2::VectorType phi = C_.Log();

    VectorType xi;
    angular_part(xi) = phi;

    double a, b;
    if (std::abs(phi) > EPSILON) {
        a = std::sin(phi) / phi;
        b = (1.0 - std::cos(phi)) / phi;
    } else  // use Taylor-series expansion
    {
        a = 1.0 - std::pow(phi, 2) / 6.0 + std::pow(phi, 4) / 120.0;
        b = phi / 2.0 - std::pow(phi, 3) / 24.0;
    }
    const Eigen::Matrix2d Vinv =
        1.0 / (a * a + b * b) * (a * Eigen::Matrix2d::Identity() - b * SO2::hat(1.0));
    linear_part(xi) = Vinv * r_;

    return xi;
}

SE2::MappingType SE2::Ad() const {
    MappingType Adj = MappingType::Identity();
    Adj.topLeftCorner<2, 2>() = C_.matrix();
    Adj(0, 2) = r_(1);
    Adj(1, 2) = -r_(0);
    return Adj;
}

SE2::MappingType SE2::Jl(const VectorType &xi) {
    const SO2::VectorType phi = angular_part(xi);
    double alpha1, alpha2;  // coefficients for ξ⋏ and (ξ⋏)²

    // use exact coefficients if |φ|<ε
    if (std::abs(phi) >
        EPSILON)  //!< @todo EPSILON may be too big since we have φ³ in the denominator
    {
        alpha1 = (1.0 - std::cos(phi)) / std::pow(phi, 2);
        alpha2 = (phi - std::sin(phi)) / std::pow(phi, 3);
    } else  // otherwise use Taylor-series expansion
    {
        double phi_2 = phi * phi;
        double phi_4 = phi_2 * phi_2;
        double phi_6 = phi_4 * phi_2;

        alpha1 = 0.5 - phi_2 / 24.0 + phi_4 / 720.0 - phi_6 / 40320.0;
        alpha2 = 1.0 / 6.0 - phi_2 / 120.0 + phi_4 / 5040.0 - phi_6 / 362880.0;
    }

    const MappingType ad_xi = ad(xi);
    return MappingType::Identity() + alpha1 * ad_xi + alpha2 * (ad_xi * ad_xi);
}

SE2::MappingType SE2::Jr(const VectorType &xi) {
    return Jl(-xi);
}

SE2::MappingType SE2::JlInverse(const VectorType &xi) {
    const SO2::VectorType phi = angular_part(xi);
    double alpha;  // coefficient for (ξ⋏)²

    // use exact expression if |φ| > ε, and φ is not an integer multiple of 2π
    if (std::abs(phi) > EPSILON &&
        std::abs(phi / (2 * M_PI) - std::round(phi / (2 * M_PI))) > EPSILON) {
        alpha = 1.0 / (phi * phi) - std::cos(phi / 2.0) / (2.0 * phi * std::sin(phi / 2.0));
    } else  // otherwise use Taylor-series expansion
    {
        alpha = 1.0 / 12.0 + std::pow(phi, 2) / 720.0 + std::pow(phi, 4) / 30240.0 +
                std::pow(phi, 6) / 1209600.0;
    }

    const MappingType ad_xi = ad(xi);
    return MappingType::Identity() - 0.5 * ad_xi + alpha * (ad_xi * ad_xi);
}

SE2::MappingType SE2::JrInverse(const VectorType &xi) {
    return JlInverse(-xi);
}

SE2::AlgebraType SE2::hat(const VectorType &xi) {
    AlgebraType x = AlgebraType::Zero();
    x.topLeftCorner<2, 2>() = SO2::hat(angular_part(xi));
    x.topRightCorner<2, 1>() = linear_part(xi);
    return x;
}

SE2::VectorType SE2::vee(const AlgebraType &x) {
    VectorType xi;
    linear_part(xi) = x.topRightCorner<2, 1>();
    angular_part(xi) = SO2::vee(x.topLeftCorner<2, 2>());
    return xi;
}

SE2 SE2::inverse() const {
    return SE2(C_.inverse(), -(C_.inverse() * r_));
}

SE2 SE2::operator*(const SE2 &rhs) const {
    return SE2(C_ * rhs.C_, C_ * rhs.r_ + r_);
}

SE2::DomainType SE2::operator*(const DomainType &x) const {
    return C_ * x + r_;
}

void SE2::setIdentity() {
    C_.setIdentity();
    r_.setZero();
}

bool SE2::isApprox(const SE2 &other) const {
    return C_.isApprox(other.C_) && r_.isApprox(other.r_);
}

bool SE2::isIdentity() const {
    return C_.isIdentity() && r_.isZero();
}

SE2::MatrixType SE2::matrix() const {
    MatrixType T = MatrixType::Identity();
    T.topLeftCorner<2, 2>() = C_.matrix();
    T.topRightCorner<2, 1>() = r_;
    return T;
}

SE2::MappingType SE2::ad(const VectorType &xi) {
    MappingType ad_xi = MappingType::Zero();

    ad_xi.topLeftCorner<2, 2>() = SO2::hat(angular_part(xi));
    ad_xi(0, 2) = xi(1);
    ad_xi(1, 2) = -xi(0);

    return ad_xi;
}

}  // namespace mange
