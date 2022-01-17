#pragma once

#include <eigen3/Eigen/Core>

namespace mange {

class SO2 {
   public:
    static constexpr int DOF = 1;
    static constexpr int DIM = 2;

    using VectorType = double;
    using MappingType = double;
    using AlgebraType = Eigen::Matrix2d;
    using MatrixType = Eigen::Matrix2d;
    using DomainType = Eigen::Matrix<double, DIM, 1>;

    SO2();
    SO2(VectorType phi);

    static SO2 Identity();
    static SO2 Random();

    static SO2 Exp(VectorType phi);
    static VectorType Log(const SO2 &X);
    static MappingType Ad(const SO2 &X);
    static MappingType Jl(VectorType phi);
    static MappingType Jr(VectorType phi);
    static MappingType JlInverse(VectorType phi);
    static MappingType JrInverse(VectorType phi);

    static AlgebraType hat(VectorType phi);
    static VectorType vee(const AlgebraType &x);

    VectorType Log() const { return Log(*this); }
    MappingType Ad() const { return Ad(*this); }

    SO2 inverse() const;

    SO2 operator*(const SO2 &rhs) const;
    DomainType operator*(const DomainType &x) const;

    void setIdentity();
    //! @todo void normalize();

    bool isApprox(const SO2 &other) const;
    bool isIdentity() const;

    const MatrixType &matrix() const { return C_; }

    // SO(2) specific methods
    double rotation() const;

   private:
    SO2(const MatrixType &C);

    MatrixType C_;
};

}  // namespace mange
