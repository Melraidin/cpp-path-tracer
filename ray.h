#ifndef RAY_H
#define RAY_H

#include <Eigen/Dense>
 
using Eigen::MatrixXd;

class Ray {
 public:
  Eigen::Vector3f o;
  Eigen::Vector3f d;

 Ray(Eigen::Vector3f o_, Eigen::Vector3f d_):
  o(o_), d(d_) {}
};

#endif
