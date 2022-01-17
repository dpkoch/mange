#pragma once

#include <eigen3/Eigen/Core>

#include "mange/SO3.hh"

namespace mange {

class SE3 {
   public:
    static constexpr int DOF = 6;
    static constexpr int DIM = 3;

    using VectorType = Eigen::Matrix<double, DOF, 1>;
    using MappingType = Eigen::Matrix<double, DOF, DOF>;
    using AlgebraType = Eigen::Matrix4d;
    using MatrixType = Eigen::Matrix4d;
    using DomainType = Eigen::Matrix<double, DIM, 1>;

    SE3();
    SE3(const VectorType &xi);
    SE3(const SO3 &rotation);
    SE3(const DomainType &translation);
    SE3(const SO3 &rotation, const DomainType &translation);

    static SE3 Identity();
    static SE3 Random();

    template <typename Derived>
    static Eigen::VectorBlock<Derived, SO3::DOF> angular_part(Eigen::MatrixBase<Derived> &vector) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, DOF);
        return vector.template tail<SO3::DOF>();
    }
    template <typename Derived>
    static const Eigen::VectorBlock<const Derived, SO3::DOF> angular_part(
        const Eigen::MatrixBase<Derived> &vector) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, DOF);
        return vector.template tail<SO3::DOF>();
    }

    template <typename Derived>
    static Eigen::VectorBlock<Derived, DOF - SO3::DOF> linear_part(
        Eigen::MatrixBase<Derived> &vector) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, DOF);
        return vector.template head<DOF - SO3::DOF>();
    }
    template <typename Derived>
    static const Eigen::VectorBlock<const Derived, DOF - SO3::DOF> linear_part(
        const Eigen::MatrixBase<Derived> &vector) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, DOF);
        return vector.template head<DOF - SO3::DOF>();
    }

    static SE3 Exp(const VectorType &xi);
    static VectorType Log(const SE3 &X);
    static MappingType Ad(const SE3 &X);
    static MappingType Jl(const VectorType &xi);
    static MappingType Jr(const VectorType &xi);
    static MappingType JlInverse(const VectorType &xi);
    static MappingType JrInverse(const VectorType &xi);

    static AlgebraType hat(const VectorType &xi);
    static VectorType vee(const AlgebraType &x);

    VectorType Log() const { return Log(*this); }
    MappingType Ad() const { return Ad(*this); }

    SE3 inverse() const;

    SE3 operator*(const SE3 &rhs) const;
    DomainType operator*(const DomainType &x) const;

    void setIdentity();
    //! @todo void normalize();

    bool isApprox(const SE3 &other) const;
    bool isIdentity() const;

    MatrixType matrix() const;

    // SE(3) specific methods
    const SO3 &C() const { return C_; }
    const DomainType &r() const { return r_; }

   private:
    static constexpr double EPSILON = 1e-12;

    static MappingType ad(const VectorType &xi);
    static SO3::MappingType jacobianQMatrix(const VectorType &xi);

    SO3 C_;
    DomainType r_;
};

}  // namespace mange
