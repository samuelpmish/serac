// Copyright (c) 2019-2020, Lawrence Livermore National Security, LLC and

// other Serac Project Developers. See the top-level LICENSE file for
// details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <array>
#include <random>
#include <vector>

struct Vector : public std::vector<double> {
  using std::vector<double>::vector;

  auto  operator()(int i) const { return this->operator[](i); }
  auto& operator()(int i) { return this->operator[](i); }

  static auto Random(int size, std::array<double, 2> range = {-1.0, 1.0})
  {
    std::default_random_engine             generator;
    std::uniform_real_distribution<double> distribution(range[0], range[1]);

    Vector x;
    x.resize(size);
    for (int i = 0; i < size; i++) {
      x(i) = distribution(generator);
    }
    return x;
  }
};

struct Matrix {
  int                 nrows_, ncols_;
  std::vector<double> values_;

  Matrix() { resize(0, 0); }
  Matrix(int rows, int cols = -1) { resize(rows, (cols == -1) ? rows : cols); }

  void resize(int rows, int cols)
  {
    nrows_ = rows;
    ncols_ = cols;
    values_.resize(rows, cols);
  }

  auto  operator()(int i, int j) const { return values_[i * ncols_ + j]; }
  auto& operator()(int i, int j) { return values_[i * ncols_ + j]; }

  static auto Random(std::array<int, 2> dimensions, std::array<double, 2> range = {-1.0, 1.0})
  {
    std::default_random_engine             generator;
    std::uniform_real_distribution<double> distribution(range[0], range[1]);

    Matrix A;
    A.resize(dimensions[0], dimensions[1]);
    for (int i = 0; i < dimensions[0]; i++) {
      for (int j = 0; j < dimensions[1]; j++) {
        A(i, j) = distribution(generator);
      }
    }
    return A;
  }
};

Matrix sym(Matrix A)
{
  Matrix Asym(A.nrows_, A.ncols_);
  for (int i = 0; i < A.nrows_; i++) {
    for (int j = 0; j < A.ncols_; j++) {
      Asym(i, j) = 0.5 * (A(i, j) + A(j, i));
    }
  }
  return Asym;
}

struct Unknown {
  struct Slice {
    Unknown&         unknowns_;
    std::vector<int> ids_;

    Slice(Unknown& unknowns, std::vector<int> ids) : unknowns_(unknowns), ids_(ids) {}

    Slice operator[](int i) { return Slice{unknowns_, {ids_[i]}}; }

    Slice operator[](std::vector<int> i)
    {
      std::vector<int> ids(i.size());
      for (size_t j = 0; j < ids.size(); j++) {
        ids[j] = ids_[i[j]];
      }
      return Slice{unknowns_, ids};
    }
  };

  struct Matvec {
    Matrix&  A_;
    Unknown& x_;
  };

  Slice operator[](int i) { return Slice{*this, {i}}; }
  Slice operator[](std::vector<int> i) { return Slice{*this, i}; }

  int size_;
};

// C * x {==, <=, >=} g
struct Constraint {
  enum class Type
  {
    Equality,
    LessThanOrEqual,
    GreaterThanOrEqual
  };

  const Matrix& C_;
  const Vector& g_;
  Type          type_;
};

auto operator==(Unknown::Matvec Ax, Vector& b) { return Constraint{Ax.A_, b, Constraint::Type::Equality}; }
auto operator<=(Unknown::Matvec Ax, Vector& b) { return Constraint{Ax.A_, b, Constraint::Type::LessThanOrEqual}; }
auto operator>=(Unknown::Matvec Ax, Vector& b) { return Constraint{Ax.A_, b, Constraint::Type::GreaterThanOrEqual}; }

auto operator==(Unknown::Slice x, const Vector& v)
{
  Matrix f;
  return Constraint{f, v, Constraint::Type::Equality};
}

auto operator*(Matrix& A, Unknown& x) { return Unknown::Matvec{A, x}; }

struct LinearSystem {
  Constraint              primary_constraint_;
  std::vector<Constraint> auxiliary_constraints_;

  LinearSystem(Constraint p) : primary_constraint_(p), auxiliary_constraints_{} {}

  void append_constraint(Constraint c) { auxiliary_constraints_.push_back(c); }
};

int main()
{
  /*
    desired interface

    // description of the physics
    Matrix K = ... ;
    Vector f = ... ;

    // constraints to be satisfied
    Matrix C = ... ;
    Vector g = ... ;

    Unknown u;

    auto sys = LinearSystem(K * u == f);
    sys.append_constraint({transpose(C) * u == g});

    auto soln = solve(sys, u);
  */

  int n = 100;
  int m = 10;

  Matrix K = sym(Matrix::Random({n, n}));
  Vector f = Vector::Random(n);

  Matrix C = Matrix::Random({m, n});
  Vector g = Vector::Random(m);

  Unknown u{n};

  auto sys = LinearSystem(K * u == f);
  sys.append_constraint(C * u == g);
  sys.append_constraint(u[{15, 32, 79, 91}] == Vector{1.0, 2.0, 3.0, 4.0});

}