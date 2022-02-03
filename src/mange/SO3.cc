#include "mange/SO3.hh"

#include <cmath>

namespace mange {

SO3::SO3() : C_(MatrixType::Identity()) {}

SO3::SO3(VectorType phi) {
    *this = Exp(phi);
}

SO3::SO3(const MatrixType &C) : C_(C) {}

SO3 SO3::Identity() {
    return SO3();
}

SO3 SO3::Random() {
    return Exp(VectorType::Random());
}

SO3 SO3::Exp(VectorType phi) {
    const double angle = phi.norm();
    const VectorType axis = phi.normalized();
    const AlgebraType axis_hat = hat(axis);

    //! @todo Replace with more standard form?
    return SO3(MatrixType(MatrixType::Identity() + std::sin(angle) * axis_hat +
                          (1.0 - std::cos(angle)) * axis_hat * axis_hat));
}

SO3::VectorType SO3::Log() const {
    if (isIdentity()) {
        return VectorType::Zero();
    }

    const double scalar = (matrix().trace() - 1.0) / 2.0;
    //! @todo This isn't stable very close to identity
    return std::acos(scalar) / (2.0 * std::sqrt(1.0 - scalar * scalar)) *
           vee((matrix() - matrix().transpose()));
}

SO3::MappingType SO3::Ad() const {
    return C_;
}

SO3::MappingType SO3::Jl(VectorType phi) {
    const double angle = phi.norm();
    const VectorType axis = phi.normalized();

    double a, b, c;
    if (angle > EPSILON) {
        const double sin_phi = std::sin(angle);

        a = sin_phi / angle;
        b = 1 - a;
        c = (1 - std::cos(angle)) / angle;
    } else {
        b = std::pow(angle, 2) / 6 - std::pow(angle, 4) / 120;
        a = 1 - b;
        c = -1 / angle + 1 + angle / 2 - std::pow(angle, 3) / 24;
    }

    return a * MappingType::Identity() + b * axis * axis.transpose() + c * hat(axis);
}

SO3::MappingType SO3::Jr(VectorType phi) {
    return Jl(-phi);
}

SO3::MappingType SO3::JlInverse(VectorType phi) {
    const double angle = phi.norm();
    const VectorType axis = phi.normalized();

    double a, b, c;
    // use exact expression if |φ| > ε, and φ is not an integer multiple of 2π
    if (angle > EPSILON &&
        std::abs(angle / (2 * M_PI) - std::round(angle / (2 * M_PI))) > EPSILON) {
        a = angle / (2 * std::tan(angle / 2));
    } else {
        a = 1 - std::pow(angle, 2) / 12 - std::pow(angle, 4) / 720;
    }
    b = 1 - a;
    c = -angle / 2;

    return a * MappingType::Identity() + b * axis * axis.transpose() + c * hat(axis);
}

SO3::MappingType SO3::JrInverse(VectorType phi) {
    return JlInverse(-phi);
}

SO3::AlgebraType SO3::hat(VectorType phi) {
    // clang-format off
    return (AlgebraType() << 0.0,   -phi(2), phi(1),
                             phi(2), 0.0,   -phi(0),
                            -phi(1), phi(0), 0.0).finished();
    // clang-format on
}

SO3::VectorType SO3::vee(const AlgebraType &x) {
    return (VectorType() << x(2, 1), x(0, 2), x(1, 0)).finished();
}

SO3 SO3::inverse() const {
    return SO3(MatrixType(C_.transpose()));
}

SO3 SO3::operator*(const SO3 &rhs) const {
    return SO3(MatrixType(C_ * rhs.C_));
}

SO3::DomainType SO3::operator*(const DomainType &x) const {
    return C_ * x;
}

void SO3::setIdentity() {
    C_.setIdentity();
}

bool SO3::isApprox(const SO3 &other) const {
    return C_.isApprox(other.C_);
}

bool SO3::isIdentity() const {
    return C_.isIdentity();
}

}  // namespace mange
