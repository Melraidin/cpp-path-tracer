#include "sphere.h"
#include "ray.h"

double Sphere::intersect(const Ray& r) const {
  auto op = this->position - r.o;
  double t = 0;
  double eps = 1e-4;
  double b = op.dot(r.d);

  // TODO Consider using Eigen's determinant().
  double det = b * b - op.dot(op) + this->radius * this->radius;
  if (det < 0) {
    return 0;
  }

  det = sqrt(det);
  return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
}
