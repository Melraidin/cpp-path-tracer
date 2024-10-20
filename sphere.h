#ifndef SPHERE_H
#define SPHERE_H

#include <Eigen/Dense>

#include "enums.h"
#include "ray.h"
 
using Eigen::MatrixXd;

class Sphere {
 public:
  double radius;
  Eigen::Vector3f position;
  Eigen::Vector3f emission;
  Eigen::Vector3f color;
  refl_t refl;

 Sphere(double r_, Eigen::Vector3f p_, Eigen::Vector3f e_, Eigen::Vector3f c_, refl_t refl_):
  radius(r_), position(p_), emission(e_), color(c_), refl(refl_) {}

  double intersect(const Ray& r) const;
};

#endif
