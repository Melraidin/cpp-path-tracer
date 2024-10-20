#include <ctime>
#include <Eigen/Dense>
#include <iostream>
// #include <math>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <png++/png.hpp>

#include "sphere.h"
#include "ray.h"

using namespace std;

Eigen::Vector3f radiance(const Ray& r, int depth, unsigned short* Xi, int E=1);
bool intersect(const Ray& r, double &t, int &id);

double uniform_random() {
  static mt19937 generator(time(nullptr));
  static uniform_real_distribution<double> distribution(0.0, 1.0);
  
  return distribution(generator);
}

Sphere spheres[] = { // Scene: radius, position, emission, color, material 
  Sphere(1e5, Eigen::Vector3f( 1e5+1,40.8,81.6), Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(.75,.25,.25), DIFF),// Left
  Sphere(1e5, Eigen::Vector3f(-1e5+99,40.8,81.6),Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(.25,.25,.75),DIFF),// Right
  Sphere(1e5, Eigen::Vector3f(50,40.8, 1e5), Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(.75,.75,.75),DIFF),//Back 
  Sphere(1e5, Eigen::Vector3f(50,40.8,-1e5+170), Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(0, 0, 0), DIFF), // Front
  Sphere(1e5, Eigen::Vector3f(50, 1e5, 81.6), Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(.75,.75,.75),DIFF), // Bottom
  Sphere(1e5, Eigen::Vector3f(50,-1e5+81.6,81.6),Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(.75,.75,.75),DIFF), // Top
  Sphere(16.5,Eigen::Vector3f(27,16.5,47), Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(1,1,1)*.999, SPEC), // Miror
  Sphere(16.5,Eigen::Vector3f(73,16.5,78), Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(1,1,1)*.999, REFR), // Glass
  Sphere(600, Eigen::Vector3f(50,681.6-.27,81.6),Eigen::Vector3f(12,12,12), Eigen::Vector3f(0, 0, 0), DIFF) // Light
};

int numSpheres = sizeof(spheres) / sizeof(Sphere);

double clamp(double x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
  return int(pow(clamp(x), 1/2.2) * 255 + 0.5);
}

int main(int argc, char* argv[]) {
  // Image dimensions.
  // int w = 512, h = 384;
  int w = 256, h = 197;
  const double field_of_view = 0.5135;
  // const double field_of_view = 0.7135;

  // Not sure why we divide the provided number of samples by 4.
  int samps = argc == 2 ? atoi(argv[1]) / 4 : 1;

  Ray cam(Eigen::Vector3f(50, 52, 295.6), Eigen::Vector3f(0, -0.042612, -1).normalized());

  auto cx = Eigen::Vector3f(w * field_of_view / h, 0, 0);

  auto crs = cx.cross(cam.d);

  cout << "cx.cross(cam.d): " << crs(0) << ", " << crs(1) << ", " << crs(2) << endl;
  crs = crs.normalized();
  cout << "crs.norm: " << crs(0) << ", " << crs(1) << ", " << crs(2) << endl;
  crs = crs * field_of_view;
  cout << "final crs: " << crs(0) << ", " << crs(1) << ", " << crs(2) << endl;

  crs = cx.cross(cam.d).normalized() * field_of_view;
  cout << "partial cy: " << crs(0) << ", " << crs(1) << ", " << crs(2) << endl;
  
  // TODO Figure out why cy differs so much.
  Eigen::Vector3f cy = cx.cross(cam.d).normalized() * field_of_view;
  cout << "real cy: " << cy(0) << ", " << cy(1) << ", " << cy(2) << endl;

  Eigen::Vector3f r;
  Eigen::Vector3f* c = new Eigen::Vector3f[w * h];

// #pragma omp parallel for schedule(dynamic, 1) private(r)
  for (int y = 0; y < h; y++) {
  // for (int y = 0; y < 2; y++) {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.0 * y / (h - 1));

    // This is the initialization vector for erand48().
    unsigned short Xi[3] = { 0, 0, (unsigned short)(y * y * y) };

    for (unsigned short x = 0; x < w; x++) {
      for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) {
	for (int sx = 0; sx < 2; sx++, r = Eigen::Vector3f(0, 0, 0)) {
	  for (int s = 0; s < samps; s++) {
	    double r1 = 2.0 * erand48(Xi);
	    double dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
	    double r2 = 2.0 * erand48(Xi);
	    double dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);

	    Eigen::Vector3f d = cx * ( ( ( sx + 0.5 + dx) / 2.0 + x) / w - 0.5) +
	      cy * ( ( ( sy + 0.5 + dy) / 2.0 + y) / h - 0.5) + cam.d;

	    // cout << "sx: " << sx << ", sy: " << sy << ", cx: " << cx << ", cy: " << cy << endl;
	    // cout << "cam.d: " << cam.d(0) << ", " << cam.d(1) << ", " << cam.d(2) << endl;
	    // cout << "dx: " << dx << ", dy: " << dy << endl;
	    // cout << "r1: " << r1 << ", r2: " << r2 << ", d: " << d(0) << ", " << d(1) << ", " << d(2) << endl;

	    // r = r + radiance(Ray(cam.o + d * 140, d.normalized()), 0, Xi) * ( 1.0 / samps );
	    r = r + radiance(Ray(cam.o + d * 140.0, d.normalized()), 0, Xi) * ( 1.0 / samps );
	    // cout << "d: " << d << ", test: " << cam.o + d * 140 << ", r: " << r << endl;
	  }

	  // cout << "r: " << r << ", clamped: " << clamp(r(0)) << ", " << clamp(r(1)) << ", " << clamp(r(2)) << endl;

	  // cout << "Multiplied: " << Eigen::Vector3f(clamp(r(0)), clamp(r(1)), clamp(r(2))) * 0.25 << endl;

	  // Push camera rays forward so they start in the interior.
	  c[i] = c[i] + Eigen::Vector3f(clamp(r(0)), clamp(r(1)), clamp(r(2))) * 0.25;

	  // cout << "c[" << i << "]: " << c[i] << endl;
	}
      }
    }
  }

  png::image<png::rgb_pixel> image(w, h);

  Eigen::Vector3f* pixel = c;
  
  for (png::uint_32 y = 0; y < h; y++) {
  // for (png::uint_32 y = 0; y < 2; y++) {
    for (png::uint_32 x = 0; x < w; x++) {
      int i = (h - y - 1) * w + x;
      pixel = &(c[i]);
      // cout << "Pixel: " << (*c)(0) << ", " << (*c)(1) << ", " << (*c)(2) << endl;
      // cout << "Pixel[" << i << "]: " << (*pixel)(0) << ", " << (*pixel)(1) << ", " << (*pixel)(2) << endl;

      auto px = png::rgb_pixel((*pixel)(0) * 255, (*pixel)(1) * 255, (*pixel)(2) * 255);
      // cout << "pixel: " << (*c)(0) * 255 << ", " << (*c)(1) * 255 << ", " << (*c)(2) * 255 << endl;

      image[h - y - 1][x] = png::rgb_pixel((*pixel)(0) * 255, (*pixel)(1) * 255, (*pixel)(2) * 255);
    }
  }

  image.write("output.png");
}

Eigen::Vector3f radiance(const Ray& r, int depth, unsigned short* Xi, int E) {
  double t;
  int id = 0;
  // cout << "ray: " << r.o(0) << ", " << r.o(1) << ", " << r.o(2) << ", d: " << r.d(0) << ", " << r.d(1) << ", " << r.d(2) << endl;
  // exit(1);
  if (!intersect(r, t, id)) {
    return Eigen::Vector3f(0, 0, 0);
  }

  if (depth > 10) {
    // cout << "Hit max depth." << endl;
    return Eigen::Vector3f(0, 0, 0);
  }

  const Sphere& obj = spheres[id];

  // Hit point.
  Eigen::Vector3f x = r.o + r.d * t;
  // cout << "Hit sphere " << id << " at time " << t << " at point " << x << endl;

  // Normal at hit point.
  Eigen::Vector3f n = (x - obj.position).normalized();
  // Normal with corrected orientation.
  Eigen::Vector3f n1 = n.dot(r.d) < 0 ? n : n * -1.0;
  // cout << "Sphere normal at hit: " << n << ", corrected: " << n1 << endl;
  Eigen::Vector3f f = obj.color;

  double p = f(0) > f(1) && f(0) > f(2) ? f(0) : f(1) > f(2) ? f(1) : f(2);
  if (++depth > 5 || !p) {
    if (erand48(Xi) < p) {
      f = f * (1.0 / p);
    } else {
      // cout << "Russian roulette failed." << endl;
      // Why multiply by E here? Is E ever non-zero?
      // return obj.emission * E;
      return obj.emission;
    }
  }

  // TODO Remove, this is for testing the initial ray casting.
  // return obj.color;

  Eigen::Vector3f w = n1;

  if (obj.refl == DIFF) {
    cout << "Hit sphere " << id << " at time " << t << " at point " << x << endl;
    cout << "Ray: " << r.o << ", " << r.d << endl;
    cout << "Xi: " << Xi << endl;
    exit(1);
    
    double r1 = 2.0 * M_PI * erand48(Xi);
    double r2 = erand48(Xi);
    double r2s = sqrt(r2);
    
    // Ideal diffuse reflection.
    Eigen::Vector3f u = ((fabs(w(0)) > 0.1 ? Eigen::Vector3f(0, 1.0, 0) : Eigen::Vector3f(1.0, 0, 0)).cross(w)).normalized();
    Eigen::Vector3f v = w.cross(u);

    // cout << "u: " << u << ", u * 2: " << u * 2.0 << endl;
    Eigen::Vector3f d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)).normalized();
    return obj.emission + f.cwiseProduct(radiance(Ray(x, d), depth, Xi));
    // cout << "Returning diffuse: " << obj.emission + f.cwiseProduct(radiance(Ray(x, d), depth, Xi)) << endl;

    // TODO This part is only in the slideshow, sort out what exactly it does.
//     Eigen::Vector3f e(0, 0, 0);
//     for (int i = 0; i < numSpheres; i++) {
//       const Sphere& s = spheres[i];
//       if (s.emission(0) <= 0 && s.emission(1) <= 0 && s.emission(2) <= 0) {
// 	// Skip objects with no emission.
// 	continue;
//       }

//       Eigen::Vector3f sw = s.position - x;
//       Eigen::Vector3f su = ((fabs(sw(1)) > 0.1 ? Eigen::Vector3f(0, 1, 0) : Eigen::Vector3f(1, 0, 0)).cross(sw)).normalized();
//       Eigen::Vector3f sv = sw.cross(su);
// n
//       double cos_a_max = sqrt(1 - s.radius * s.radius / (x - s.position).dot(x - s.position));
//       double eps1 = erand48(Xi);
//       double eps2 = erand48(Xi);
//       double cos_a = 1 - eps1 + eps1 * cos_a_max;
//       double sin_a = sqrt(1 - cos_a * cos_a);
//       double phi = 2 * M_PI * eps2;
//       Eigen::Vector3f l = (su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a).normalized();

//       if (intersect(Ray(x, l), t, id) && id == i) {
// 	double omega = 2 * M_PI * (1 - cos_a_max);
// 	e = e + f.cwiseProduct(s.emission * l.dot(n1) * omega) * M_1_PI;
//       }
//     }
//
//     return obj.emission * E + e + f.cwiseProduct(radiance(Ray(x, d), depth, Xi, 0));
  } else if (obj.refl == SPEC) {
    // Ideal specular reflection.
    return obj.emission + f.cwiseProduct(radiance(Ray(x, r.d - n * 2.0 * n.dot(r.d)), depth, Xi));
  }

  // Refraction.
  Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // Ideal dielectric REFRACTION 
  bool into = n.dot(n1) > 0;                // Ray from outside going in? 
  double nc = 1.0;
  double nt = 1.5;
  double nnt = into ? nc / nt : nt / nc;
  double ddn = r.d.dot(n1);
  double cos2t; 
  if ((cos2t = 1.0 - nnt * nnt * ( 1.0 - ddn * ddn )) < 0) {    // Total internal reflection 
    return obj.emission + f.cwiseProduct(radiance(reflRay, depth, Xi));
  }
  Eigen::Vector3f tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1.0-(into?-ddn:tdir.dot(n)); 
  double Re=R0+(1.0-R0)*c*c*c*c*c,Tr=1.0-Re,P=0.25 + 0.5 * Re;
  double RP = Re / P;
  double TP = Tr / (1.0 - P);

  if (depth > 2) {
    Eigen::Vector3f opt1;
    if (erand48(Xi) < P) {
      opt1 = radiance(reflRay, depth, Xi) * RP;
    } else {
      opt1 = radiance(Ray(x, tdir), depth, Xi) * TP;
    }

    return obj.emission + f.cwiseProduct(opt1);
  } else {
    return obj.emission + f.cwiseProduct(radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
  }
  // return obj.emission + f.cwiseProduct(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette 
  // 		    radiance(reflRay,depth,Xi)*RP : radiance(Ray(x,tdir),depth,Xi)*TP) : 
  // 			radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
}

bool intersect(const Ray& r, double &t, int &id) {
  double d;
  double inf = t = 1e20;

  for (int i = int(numSpheres); i--;) {
    if ((d = spheres[i].intersect(r)) && d < t) {
      t = d;
      id = i;
    }
  }

  return t < inf;
}
