#ifndef MANGE_SE2_H
#define MANGE_SE2_H

#include <mange/SO2.h>

#include <eigen3/Eigen/Core>

namespace mange
{

class SE2
{
public:
  SE2();
  SE2(const Eigen::Vector3d &xi);

  static SE2 Identity();
  static SE2 Random();

  static SE2 Exp(const Eigen::Vector3d &xi);
  static Eigen::Vector3d Log(const SE2 &X);
  static Eigen::Matrix3d Ad(const SE2 &X);
  static Eigen::Matrix3d Jl(const Eigen::Vector3d &xi);
  static Eigen::Matrix3d Jr(const Eigen::Vector3d &xi);
  static Eigen::Matrix3d JlInverse(const Eigen::Vector3d &xi);
  static Eigen::Matrix3d JrInverse(const Eigen::Vector3d &xi);

  static Eigen::Matrix3d hat(const Eigen::Vector3d &xi);
  static Eigen::Vector3d vee(const Eigen::Matrix3d &x);

  Eigen::Vector3d Log() const { return Log(*this); }
  Eigen::Matrix3d Ad() const { return Ad(*this); }

  SE2 inverse() const;

  SE2 operator*(const SE2 &rhs) const;
  Eigen::Vector2d operator*(const Eigen::Vector2d &x) const;

  void setIdentity();
  //! @todo void normalize();

  bool isApprox(const SE2 &other) const;
  bool isIdentity() const;

  Eigen::Matrix3d matrix() const;

  // SE(2) specific methods
  const SO2 &C() const { return C_; }
  const Eigen::Vector2d &r() const { return r_; }

private:
  static constexpr double EPSILON = 1e-12;

  SE2(const SO2 &C, const Eigen::Vector2d &r);
  static Eigen::Matrix3d ad(const Eigen::Vector3d &xi);

  SO2 C_;
  Eigen::Vector2d r_;
};

} // namespace mange

#endif // MANGE_SE2_H
