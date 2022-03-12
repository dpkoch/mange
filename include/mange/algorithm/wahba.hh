#pragma once

#include <vector>

#include <eigen3/Eigen/Dense>

namespace mange {
class SO3;

namespace algorithm {

struct WahbaSample {
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double weight = 1.0;
};

SO3 wahba(const std::vector<WahbaSample>& samples);

SO3 wahba(const Eigen::Vector3d& w, const Eigen::Vector3d& v);

}  // namespace algorithm
}  // namespace mange
