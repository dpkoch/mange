#include <mange/SE2.h>

namespace mange
{

SE2::SE2()
{
  //! @todo
}

SE2::SE2(const Eigen::Vector3d &phi)
{
  //! @todo
}

SE2 SE2::Identity()
{
  //! @todo
  return SE2();
}

SE2 SE2::Random()
{
  //! @todo
  return SE2();
}

SE2 SE2::Exp(const Eigen::Vector3d &phi)
{
  //! @todo
  return SE2();
}

Eigen::Vector3d SE2::Log(const SE2 &X)
{
  //! @todo
  return Eigen::Vector3d();
}

Eigen::Matrix3d SE2::Ad(const SE2 &X)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::Jl(const Eigen::Vector3d &phi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::Jr(const Eigen::Vector3d &phi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::JlInverse(const Eigen::Vector3d &phi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::JrInverse(const Eigen::Vector3d &phi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Matrix3d SE2::hat(const Eigen::Vector3d &phi)
{
  //! @todo
  return Eigen::Matrix3d();
}

Eigen::Vector3d SE2::vee(const Eigen::Matrix3d &x)
{
  //! @todo
}

SE2 SE2::inverse() const
{
  //! @todo
  return SE2();
}

SE2 SE2::operator*(const SE2 &rhs) const
{
  //! @todo
  return SE2();
}

Eigen::Vector2d SE2::operator*(const Eigen::Vector2d &x) const
{
  //! @todo
  return Eigen::Vector2d();
}

void SE2::setIdentity()
{
  //! @todo
}

bool SE2::isApprox(const SE2 &other) const
{
  //! @todo
  return false;
}

bool SE2::isIdentity() const
{
  //! @todo
  return false;
}

Eigen::Matrix3d SE2::matrix() const
{
  //! @todo
  return Eigen::Matrix3d();
}

} // namespace mange
