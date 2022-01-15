#ifndef MANGE_SE2_H
#define MANGE_SE2_H

#include <eigen3/Eigen/Core>

#include "mange/SO2.h"

namespace mange {

class SE2 {
   public:
    static constexpr int DOF = 3;
    static constexpr int DIM = 2;

    using VectorType = Eigen::Matrix<double, DOF, 1>;
    using MappingType = Eigen::Matrix<double, DOF, DOF>;
    using AlgebraType = Eigen::Matrix<double, DOF, DOF>;
    using MatrixType = Eigen::Matrix3d;
    using DomainType = Eigen::Matrix<double, DIM, 1>;

    SE2();
    SE2(const VectorType &xi);
    SE2(const SO2 &rotation);
    SE2(const DomainType &translation);
    SE2(const SO2 &rotation, const DomainType &translation);

    static SE2 Identity();
    static SE2 Random();

    static SE2 Exp(const VectorType &xi);
    static VectorType Log(const SE2 &X);
    static MappingType Ad(const SE2 &X);
    static MappingType Jl(const VectorType &xi);
    static MappingType Jr(const VectorType &xi);
    static MappingType JlInverse(const VectorType &xi);
    static MappingType JrInverse(const VectorType &xi);

    static AlgebraType hat(const VectorType &xi);
    static VectorType vee(const AlgebraType &x);

    VectorType Log() const { return Log(*this); }
    MappingType Ad() const { return Ad(*this); }

    SE2 inverse() const;

    SE2 operator*(const SE2 &rhs) const;
    DomainType operator*(const DomainType &x) const;

    void setIdentity();
    //! @todo void normalize();

    bool isApprox(const SE2 &other) const;
    bool isIdentity() const;

    MatrixType matrix() const;

    // SE(2) specific methods
    const SO2 &C() const { return C_; }
    const Eigen::Vector2d &r() const { return r_; }

   private:
    static constexpr double EPSILON = 1e-12;

    static Eigen::Matrix3d ad(const Eigen::Vector3d &xi);

    SO2 C_;
    Eigen::Vector2d r_;
};

}  // namespace mange

#endif  // MANGE_SE2_H
