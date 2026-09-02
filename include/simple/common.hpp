#pragma once

#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Dense>

namespace simple {

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using MatX = Eigen::MatrixXd;

using Vec3Arr = std::vector<Vec3>;
using Vec2Arr = std::vector<Vec2>;
using VecDArr = std::vector<double>;

constexpr double PI = 3.14159265358979323846;
constexpr double EPS = 1e-12;

enum class ParamDir { U = 0, V = 1 };

inline double degToRad(double deg) { return deg * PI / 180.0; }
inline double radToDeg(double rad) { return rad * 180.0 / PI; }

inline double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

struct IterVersion {
    int version;
    std::string label;
    std::vector<Vec3Arr> ctrlCurvesC0;
    std::vector<Vec3Arr> ctrlCurvesC1;
};

} // namespace simple
