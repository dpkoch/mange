#ifndef MANGE_SO2_H
#define MANGE_SO2_H

#include <eigen3/Eigen/Core>

namespace mange
{

class SO2
{
public:
  SO2();
  SO2(double rotation);

  static SO2 Identity();
  static SO2 Random();

  static SO2 Exp(double phi);
  static double Log(const SO2 &X);
  static double Ad(const SO2 &X);
  static double Jl(double phi);
  static double Jr(double phi);
  static double JlInverse(double phi);
  static double JrInverse(double phi);

  static Eigen::Matrix2d hat(double phi);
  static double vee(const Eigen::Matrix2d &x);

  double Log() const { return Log(*this); }
  double Ad() const { return Ad(*this); }

  SO2 inverse() const;

  SO2 operator*(const SO2 &rhs) const;
  Eigen::Vector2d operator*(const Eigen::Vector2d &x) const;

  void setIdentity();
  //! @todo void normalize();

  bool isApprox(const SO2 &other) const;
  bool isIdentity() const;

  const Eigen::Matrix2d &matrix() const { return C_; }

  // SO(2) specific methods
  double rotation() const;

private:
  SO2(const Eigen::Matrix2d &C);

  Eigen::Matrix2d C_;
};

} // namespace mange

#endif // MANGE_SO2_H
