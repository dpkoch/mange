#ifndef MANGE_SO3_H
#define MANGE_SO3_H

#include <tuple>

#include <eigen3/Eigen/Core>

namespace mange {

class SO3 {
   public:
    static constexpr int DOF = 3;
    static constexpr int DIM = 3;

    using VectorType = Eigen::Matrix<double, DOF, 1>;
    using MappingType = Eigen::Matrix<double, DOF, DOF>;
    using AlgebraType = Eigen::Matrix3d;
    using MatrixType = Eigen::Matrix3d;
    using DomainType = Eigen::Matrix<double, DIM, 1>;

    SO3();
    SO3(VectorType phi);

    static SO3 Identity();
    static SO3 Random();

    static SO3 Exp(VectorType phi);
    static VectorType Log(const SO3 &X);
    static MappingType Ad(const SO3 &X);
    static MappingType Jl(VectorType phi);
    static MappingType Jr(VectorType phi);
    static MappingType JlInverse(VectorType phi);
    static MappingType JrInverse(VectorType phi);

    static AlgebraType hat(VectorType phi);
    static VectorType vee(const AlgebraType &x);

    VectorType Log() const { return Log(*this); }
    MappingType Ad() const { return Ad(*this); }

    SO3 inverse() const;

    SO3 operator*(const SO3 &rhs) const;
    DomainType operator*(const DomainType &x) const;

    void setIdentity();
    //! @todo void normalize();

    bool isApprox(const SO3 &other) const;
    bool isIdentity() const;

    const MatrixType &matrix() const { return C_; }

   private:
    static constexpr double EPSILON = 1e-12;

    SO3(const MatrixType &C);

    MatrixType C_;
};

}  // namespace mange

#endif  // MANGE_SO3_H
