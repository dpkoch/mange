#include "mange/algorithm/wahba.hh"

#include <eigen3/Eigen/SVD>

#include "mange/SO3.hh"

namespace mange::algorithm {

SO3 wahba(const std::vector<WahbaSample>& samples) {
    // Algorithm from https://en.wikipedia.org/wiki/Wahba%27s_problem
    using Matrix = Eigen::Matrix<double, SO3::DIM, SO3::DIM>;

    Matrix B = Matrix::Zero();
    for (const auto& sample : samples) {
        B += sample.weight * sample.w * sample.v.transpose();
    }

    Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Matrix& U = svd.matrixU();
    const Matrix& V = svd.matrixV();

    const Matrix M =
        (SO3::DomainType() << 1.0, 1.0, U.determinant() * V.determinant()).finished().asDiagonal();

    return SO3(U * M * V.transpose());
}

SO3 wahba(const Eigen::Vector3d& w, const Eigen::Vector3d& v) {
    return wahba(std::vector<WahbaSample>{{w, v}});
}

}  // namespace mange::algorithm