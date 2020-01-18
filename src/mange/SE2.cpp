#include <mange/SE2.h>

namespace mange
{

SE2::SE2() :
  r_(Eigen::Vector2d::Zero())
{
}

SE2::SE2(const Eigen::Vector3d &xi)
{
  *this = Exp(xi);
}

SE2::SE2(const SO2 &C, const Eigen::Vector2d &r) :
  C_(C),
  r_(r)
{
}

SE2 SE2::Identity()
{
  return SE2();
}

SE2 SE2::Random()
{
  return SE2(SO2::Random(), Eigen::Vector2d::Random());
}

SE2 SE2::Exp(const Eigen::Vector3d &xi)
{
  SE2 result;

  double phi = xi(2);
  result.C_ = SO2(phi);

  double a, b;
  if (std::abs(phi) > EPSILON)
  {
    a = std::sin(phi) / phi;
    b = (1.0 - std::cos(phi)) / phi;
  }
  else // use Taylor-series expansion
  {
    a = 1.0 - std::pow(phi, 2)/6.0 + std::pow(phi, 4)/120.0;
    b = phi/2.0 - std::pow(phi, 3)/24.0;
  }
  Eigen::Matrix2d V = (a * Eigen::Matrix2d::Identity() + b * SO2::hat(1.0));
  result.r_ =  V * xi.topRows<2>();

  return result;
}

Eigen::Vector3d SE2::Log(const SE2 &X)
{
  double phi = X.C_.Log();

  Eigen::Vector3d xi;
  xi(2) = phi;

  double a, b;
  if (std::abs(phi) > EPSILON)
  {
    a = std::sin(phi) / phi;
    b = (1.0 - std::cos(phi)) / phi;
  }
  else
  {
    a = 1.0 - std::pow(phi, 2)/6.0 + std::pow(phi, 4)/120.0;
    b = phi/2.0 - std::pow(phi, 3)/24.0;
  }
  Eigen::Matrix2d Vinv = 1.0 / (a*a + b*b) * (a * Eigen::Matrix2d::Identity() - b * SO2::hat(1.0));
  xi.topRows<2>() = Vinv * X.r_;

  return xi;
}

Eigen::Matrix3d SE2::Ad(const SE2 &X)
{
  Eigen::Matrix3d Adj = Eigen::Matrix3d::Identity();
  Adj.topLeftCorner<2,2>() = X.C_.matrix();
  Adj(0,2) = X.r_(1);
  Adj(1,2) = -X.r_(0);
  return Adj;
}

Eigen::Matrix3d SE2::Jl(const Eigen::Vector3d &xi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::Jr(const Eigen::Vector3d &xi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::JlInverse(const Eigen::Vector3d &xi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::JrInverse(const Eigen::Vector3d &xi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::hat(const Eigen::Vector3d &xi)
{
  Eigen::Matrix3d x = Eigen::Matrix3d::Zero();
  x.topLeftCorner<2,2>() = SO2::hat(xi(2));
  x.topRightCorner<2,1>() = xi.topRows<2>();
  return x;
}

Eigen::Vector3d SE2::vee(const Eigen::Matrix3d &x)
{
  Eigen::Vector3d xi;
  xi.topRows<2>() = x.topRightCorner<2,1>();
  xi(2) = SO2::vee(x.topLeftCorner<2,2>());
  return xi;
}

SE2 SE2::inverse() const
{
  return SE2(C_.inverse(), -(C_.inverse()*r_));
}

SE2 SE2::operator*(const SE2 &rhs) const
{
  return SE2(C_*rhs.C_, C_*rhs.r_ + r_);
}

Eigen::Vector2d SE2::operator*(const Eigen::Vector2d &x) const
{
  return C_*x + r_;
}

void SE2::setIdentity()
{
  C_.setIdentity();
  r_.setZero();
}

bool SE2::isApprox(const SE2 &other) const
{
  return C_.isApprox(other.C_) && r_.isApprox(other.r_);
}

bool SE2::isIdentity() const
{
  return C_.isIdentity() && r_.isZero();
}

Eigen::Matrix3d SE2::matrix() const
{
  Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
  T.topLeftCorner<2,2>() = C_.matrix();
  T.topRightCorner<2,1>() = r_;
  return T;
}

} // namespace mange
