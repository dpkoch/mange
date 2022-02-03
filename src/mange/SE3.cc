#include "mange/SE3.hh"

namespace mange {

SE3::SE3() : r_(DomainType::Zero()) {}

SE3::SE3(const VectorType &xi) {
    *this = Exp(xi);
}

SE3::SE3(const SO3 &rotation) : C_(rotation) {}

SE3::SE3(const DomainType &translation) : r_(translation) {}

SE3::SE3(const SO3 &rotation, const DomainType &translation) : C_(rotation), r_(translation) {}

SE3 SE3::Identity() {
    return SE3();
}

SE3 SE3::Random() {
    return SE3(SO3::Random(), DomainType::Random());
}

SE3 SE3::Exp(const VectorType &xi) {
    SE3 result;
    result.C_ = SO3::Exp(angular_part(xi));
    result.r_ = SO3::Jl(angular_part(xi)) * linear_part(xi);
    return result;
}

SE3::VectorType SE3::Log() const {
    VectorType result;
    angular_part(result) = C_.Log();
    linear_part(result) = SO3::JlInverse(angular_part(result)) * r_;
    return result;
}

SE3::MappingType SE3::Ad() const {
    MappingType result = MappingType::Zero();
    result.topLeftCorner<3, 3>() = C_.matrix();
    result.bottomRightCorner<3, 3>() = C_.matrix();
    result.topRightCorner<3, 3>() = SO3::hat(r_) * C_.matrix();
    return result;
}

SE3::MappingType SE3::Jl(const VectorType &xi) {
    const SO3::MappingType J = SO3::Jl(angular_part(xi));
    const SO3::MappingType Q = jacobianQMatrix(xi);

    MappingType result = MappingType::Zero();
    result.topLeftCorner<SO3::DOF, SO3::DOF>() = J;
    result.bottomRightCorner<SO3::DOF, SO3::DOF>() = J;
    result.topRightCorner<3, 3>() = Q;
    return result;
}

SE3::MappingType SE3::Jr(const VectorType &xi) {
    return Jl(-xi);
}

SE3::MappingType SE3::JlInverse(const VectorType &xi) {
    const SO3::MappingType Jinverse = SO3::JlInverse(angular_part(xi));
    const SO3::MappingType Q = jacobianQMatrix(xi);

    MappingType result = MappingType::Zero();
    result.topLeftCorner<SO3::DOF, SO3::DOF>() = Jinverse;
    result.bottomRightCorner<SO3::DOF, SO3::DOF>() = Jinverse;
    result.topRightCorner<3, 3>() = -Jinverse * Q * Jinverse;
    return result;
}

SE3::MappingType SE3::JrInverse(const VectorType &xi) {
    return JlInverse(-xi);
}

SE3::AlgebraType SE3::hat(const VectorType &xi) {
    AlgebraType x = AlgebraType::Zero();
    x.topLeftCorner<DIM, DIM>() = SO3::hat(angular_part(xi));
    x.topRightCorner<DIM, 1>() = linear_part(xi);
    return x;
}

SE3::VectorType SE3::vee(const AlgebraType &x) {
    VectorType xi;
    angular_part(xi) = SO3::vee(x.topLeftCorner<DIM, DIM>());
    linear_part(xi) = x.topRightCorner<DIM, 1>();
    return xi;
}

SE3 SE3::inverse() const {
    return SE3(C_.inverse(), -(C_.inverse() * r_));
}

SE3 SE3::operator*(const SE3 &rhs) const {
    return SE3(C_ * rhs.C_, C_ * rhs.r_ + r_);
}

SE3::DomainType SE3::operator*(const DomainType &x) const {
    return C_ * x + r_;
}

void SE3::setIdentity() {
    C_.setIdentity();
    r_.setZero();
}

bool SE3::isApprox(const SE3 &other) const {
    return C_.isApprox(other.C_) && r_.isApprox(other.r_);
}

bool SE3::isIdentity() const {
    return C_.isIdentity() && r_.isZero();
}

SE3::MatrixType SE3::matrix() const {
    MatrixType T = MatrixType::Identity();
    T.topLeftCorner<DIM, DIM>() = C_.matrix();
    T.topRightCorner<DIM, 1>() = r_;
    return T;
}

SE3::MappingType SE3::ad(const VectorType &xi) {
    MappingType result = MappingType::Zero();
    result.topLeftCorner<SO3::DOF, SO3::DOF>() = SO3::hat(angular_part(xi));
    result.bottomRightCorner<SO3::DOF, SO3::DOF>() = result.topLeftCorner<SO3::DOF, SO3::DOF>();
    result.topRightCorner<SO3::DOF, SO3::DOF>() = SO3::hat(linear_part(xi));
    return result;
}

SO3::MappingType SE3::jacobianQMatrix(const VectorType &xi) {
    const double angle = angular_part(xi).norm();
    const SO3::VectorType axis = angular_part(xi).normalized();

    const double a = 0.5;
    double b, c, d;
    const double angle_squared = std::pow(angle, 2);
    if (std::abs(angle) > EPSILON) {
        const double sin_angle = std::sin(angle);
        const double cos_angle = std::cos(angle);

        b = (angle - sin_angle) / angle_squared;
        c = (0.5 * angle_squared + cos_angle - 1) / angle_squared;
        d = (angle - 1.5 * sin_angle + angle * 0.5 * cos_angle) / angle_squared;
    } else {
        const double angle_cubed = std::pow(angle, 3);
        b = angle / 6 - angle_cubed / 120;
        c = angle_squared / 24 - std::pow(angle, 4) / 720;
        d = angle_cubed / 120;
    }

    const SO3::AlgebraType axis_hat = SO3::hat(axis);
    const SO3::AlgebraType axis_hat_squared = axis_hat * axis_hat;
    const SO3::AlgebraType rho_hat = SO3::hat(linear_part(xi));

    return a * rho_hat +
           b * (axis_hat * rho_hat + rho_hat * axis_hat + angle * axis_hat * rho_hat * axis_hat) +
           c * (axis_hat_squared * rho_hat + rho_hat * axis_hat_squared -
                3 * axis_hat * rho_hat * axis_hat) +
           d * (axis_hat * rho_hat * axis_hat_squared + axis_hat_squared * rho_hat * axis_hat);
}

}  // namespace mange
