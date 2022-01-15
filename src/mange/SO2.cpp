#include "mange/SO2.h"

#include <random>

namespace mange {

SO2::SO2() : C_(MatrixType::Identity()) {}

SO2::SO2(VectorType phi) {
    *this = Exp(phi);
}

SO2::SO2(const MatrixType &C) : C_(C) {}

SO2 SO2::Identity() {
    return SO2();
}

SO2 SO2::Random() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dist(-100.0, 100.0);

    return SO2(dist(gen));
}

SO2 SO2::Exp(VectorType phi) {
    SO2 result;
    result.C_ << std::cos(phi), -std::sin(phi), std::sin(phi), std::cos(phi);
    return result;
}

SO2::VectorType SO2::Log(const SO2 &X) {
    return std::atan2(X.C_(1, 0), X.C_(0, 0));
}

SO2::MappingType SO2::Ad(const SO2 &X) {
    return 1.0;
}

SO2::MappingType SO2::Jl(VectorType phi) {
    return 1.0;
}

SO2::MappingType SO2::Jr(VectorType phi) {
    return 1.0;
}

SO2::MappingType SO2::JlInverse(VectorType phi) {
    return 1.0;
}

SO2::MappingType SO2::JrInverse(VectorType phi) {
    return 1.0;
}

SO2::AlgebraType SO2::hat(VectorType phi) {
    return (AlgebraType() << 0.0, -phi, phi, 0.0).finished();
}

SO2::VectorType SO2::vee(const AlgebraType &x) {
    return x(1, 0);
}

SO2 SO2::inverse() const {
    return SO2(C_.transpose());
}

SO2 SO2::operator*(const SO2 &rhs) const {
    return SO2(C_ * rhs.C_);
}

SO2::DomainType SO2::operator*(const DomainType &x) const {
    return C_ * x;
}

double SO2::rotation() const {
    return Log();
}

void SO2::setIdentity() {
    C_.setIdentity();
}

bool SO2::isApprox(const SO2 &other) const {
    return C_.isApprox(other.C_);
}

bool SO2::isIdentity() const {
    return C_.isIdentity();
}

}  // namespace mange
