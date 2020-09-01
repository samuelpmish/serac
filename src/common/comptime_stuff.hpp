#include "mfem.hpp"

#include <array>


template <size_t rows, size_t cols>
class FixedMatrix
{
  public:
    constexpr double& operator()(size_t i, size_t j) { return data_[i + (rows * j)]; }
    double* data() {return data_; }
  private:
    std::array<double, rows * cols> data_ = {};
};

struct FixedIP
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double weight = 0.0;
};

template <int points>
constexpr std::array<FixedIP, points> gauss_legendre() {
  std::array<FixedIP, points> result;
  switch (points) {
    case 1:
     result[0].x = 0.5;
     result[0].weight = 1.0;
     break;
  }
  return result;
}

// Assuming H1 elements for now...
template <size_t dof, size_t dim, mfem::Geometry::Type geom, size_t order>
struct FixedElement
{
  constexpr FixedMatrix<dof, dim> dShape(const FixedIP& ip) {
    FixedMatrix<dof, dim> result;
    std::array<double, order + 1> dshape_dx{};
    // Barycentric
    if (order == 0) {
      dshape_dx[0] = 0.0;
    } else {
      double lk = 1.0;
      for (auto k = 0; k < order; k++) {
        // Spaghetti
        lk *= ip.x /*- x[k] */;
      }
    }
    result(0,0) = 0;
    result(1,0) = 1;
    for (auto i = 0; i < order; i++) {
      result(i+1, 0) = dshape_dx[i];
    }

    return result;
  }
};


template <mfem::Geometry::Type geom, int order>
struct FixedIntRules {};

template <int order>
struct FixedIntRules<mfem::Geometry::Type::SEGMENT, order>
{
  static constexpr int points = order / 2 + 1;
  constexpr std::array<FixedIP, points> operator()() {return gauss_legendre<points>();}
};
