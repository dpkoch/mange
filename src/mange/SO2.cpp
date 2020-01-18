#include <mange/SO2.h>

#include <random>
namespace mange
{

SO2::SO2() :
  C_(Eigen::Matrix2d::Identity())
{
}

SO2::SO2(double rotation)
{
  *this = Exp(rotation);
}

SO2::SO2(const Eigen::Matrix2d &C) :
  C_(C)
{
}

SO2 SO2::Identity()
{
  return SO2();
}

SO2 SO2::Random()
{
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<> dist(-100.0, 100.0);

  return SO2(dist(gen));
}

SO2 SO2::Exp(double phi)
{
  Eigen::Matrix2d C;
  C << std::cos(phi), -std::sin(phi),
       std::sin(phi),  std::cos(phi);
  return SO2(C);
}

double SO2::Log(const SO2 &X)
{
  return std::atan2(X.C_(1,0), X.C_(0,0));
}

double SO2::Ad(const SO2 &X)
{
  return 1.0;
}

double SO2::Jl(double phi)
{
  return 1.0;
}

double SO2::Jr(double phi)
{
  return 1.0;
}

double SO2::JlInverse(double phi)
{
  return 1.0;
}

double SO2::JrInverse(double phi)
{
  return 1.0;
}

Eigen::Matrix2d SO2::hat(double phi)
{
  Eigen::Matrix2d x;
  x << 0.0, -phi,
       phi,  0.0;
  return x;
}

double SO2::vee(const Eigen::Matrix2d &x)
{
  return x(1,0);
}

SO2 SO2::inverse() const
{
  return SO2(C_.transpose());
}

SO2 SO2::operator*(const SO2 &rhs) const
{
  return SO2(C_ * rhs.C_);
}

Eigen::Vector2d SO2::operator*(const Eigen::Vector2d &x) const
{
  return C_ * x;
}

double SO2::rotation() const
{
  return Log();
}

void SO2::setIdentity()
{
  C_.setIdentity();
}

bool SO2::isApprox(const SO2 &other) const
{
  return C_.isApprox(other.C_);
}

bool SO2::isIdentity() const
{
  return C_.isIdentity();
}

} // namespace mange
