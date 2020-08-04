// Copyright (c) 2019-2020, Lawrence Livermore National Security, LLC and

// other Serac Project Developers. See the top-level LICENSE file for
// details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <variant>
#include <vector>

#include "mfem.hpp"
#include "serac_config.hpp"

struct coo_entry {
  int    row;
  int    col;
  double value;
  bool   operator<(const coo_entry& other) const { return (row < other.row) || (row == other.row && col < other.col); }
};

// mfem::SparseMatrix doesn't have a COO constructor (?)
mfem::SparseMatrix sparse_matrix_from_tuples(std::vector<coo_entry> tuples, int nrows, int ncols)
{
  size_t n = tuples.size();

  std::sort(std::begin(tuples), std::end(tuples));

  int*    row_ptr = new int[nrows + 1];
  int*    col_ind = new int[n];
  double* values  = new double[n];

  // reduce sorted COO entries by key
  int nnz      = -1;
  int prev_row = 0;
  int prev_col = -1;
  row_ptr[0]   = 0;
  for (auto [row, col, value] : tuples) {
    if ((prev_row != row) || (prev_col != col)) {
      nnz++;
      col_ind[nnz] = col;
      values[nnz]  = value;
      for (int r = prev_row + 1; r <= row; r++) {
        row_ptr[r] = nnz;
      }
      prev_row = row;
      prev_col = col;
    } else {
      values[nnz] += value;
    }
  }
  row_ptr[nrows] = nnz;

  // mfem::SparseMatrix ctor assumes ownership,
  // so we don't deallocate rows, cols, or values
  return mfem::SparseMatrix(row_ptr, col_ind, values, nrows, ncols);
}

void must_contain(std::string haystack, std::string needle)
{
  if (haystack.find(needle) == std::string::npos) {
    std::cout << "error: \"" << haystack << "\" does not contain \"" << needle << "\"\n";
    exit(1);
  }
}

std::string find(std::string haystack, std::vector<std::string> needles)
{
  for (auto needle : needles) {
    if (haystack.find(needle) == std::string::npos) {
      return needle;
    }
  }
  return "";
}

mfem::SparseMatrix read_matrix_market(std::string filename)
{
  std::ifstream infile(filename);
  if (infile) {
    std::cout << "reading matrix from file: " << filename << std::endl;

    std::string header;
    std::getline(infile, header);

    // only allow real-valued sparse matrices
    must_contain(header, "real");
    must_contain(header, "matrix");
    must_contain(header, "coordinate");

    // detect the type of symmetry, if any
    std::string type = find(header, {"general", "symmetric", "skew-symmetric"});

    coo_entry              entry;
    std::vector<coo_entry> tuples;
    int                    nrows, ncols, nnz;

    std::string line;
    bool        uninitialized = true;
    while (std::getline(infile, line)) {
      if (line[0] == '%') {
        continue;  // disregard comments
      } else {
        std::stringstream ss(line);
        if (uninitialized) {
          ss >> nrows;
          ss >> ncols;
          ss >> nnz;
          tuples.reserve(2 * nnz);
          uninitialized = false;
        } else {
          ss >> entry.row;
          ss >> entry.col;
          ss >> entry.value;
          tuples.push_back(entry);

          if (type == "symmetric") {
            tuples.push_back({entry.col, entry.row, entry.value});
          }

          if (type == "skew-symmetric") {
            tuples.push_back({entry.col, entry.row, -entry.value});
          }
        }
      }
    }

    return sparse_matrix_from_tuples(tuples, nrows, ncols);
  } else {
    std::cout << "file: " << filename << "not found. Exiting..." << std::endl;
    exit(1);
  }
}

mfem::SparseMatrix read_sparse_matrix(std::string filename)
{
  std::ifstream infile(filename);
  if (infile) {
    std::cout << "reading matrix from file: " << filename << std::endl;

    std::vector<coo_entry> tuples;

    coo_entry entry;
    int       nrows, ncols, unused;
    infile >> unused;
    infile >> nrows;
    infile >> unused;
    infile >> ncols;

    // HypreParMatrix::Print() prints 1 less than the number of rows/cols (?)
    nrows++;
    ncols++;

    std::string line;
    while (std::getline(infile, line)) {
      std::stringstream ss(line);
      ss >> entry.row;
      ss >> entry.col;
      ss >> entry.value;
      tuples.push_back(entry);
    }

    return sparse_matrix_from_tuples(tuples, nrows, ncols);
  } else {
    std::cout << "file: " << filename << "not found. Exiting..." << std::endl;
    exit(1);
  }
}

mfem::Vector random_vector(int n)
{
  static std::default_random_engine             generator;
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  mfem::Vector                                  g(n);

  for (int i = 0; i < n; i++) {
    g[i] = distribution(generator);
  }

  return g;
}

mfem::SparseMatrix random_sparse_matrix(int nrows, int ncols, double density)
{
  // clip density value to be in the interval [~0.0, 1.0]
  double clipped_density = std::max(1.0e-8, std::min(density, 1.0));

  double                                        max_stepsize = (2.0 / clipped_density) - 1.0;
  static std::default_random_engine             generator;
  static std::uniform_real_distribution<double> distribution(0.0, 1.0);
  static std::uniform_real_distribution<double> step(1.0, max_stepsize);

  size_t nnz_estimate = nrows * ncols * clipped_density;

  std::vector<int>    row_ptr(nrows + 1);
  std::vector<int>    col_ind;
  std::vector<double> values;

  col_ind.reserve(2 * nnz_estimate);
  values.reserve(2 * nnz_estimate);

  row_ptr[0] = 0;
  for (int row = 0; row < nrows; row++) {
    double col = step(generator) - 1.0;
    while (col < ncols) {
      col_ind.push_back(floor(col));
      values.push_back(distribution(generator));
      col += step(generator);
    }
    row_ptr[row + 1] = values.size();
  }

  int nnz = row_ptr[nrows];

  int*    row_ptr_raw = new int[nrows + 1];
  int*    col_ind_raw = new int[nnz];
  double* values_raw  = new double[nnz];

  std::copy(row_ptr.begin(), row_ptr.end(), row_ptr_raw);
  std::copy(col_ind.begin(), col_ind.end(), col_ind_raw);
  std::copy(values.begin(), values.end(), values_raw);

  // mfem::SparseMatrix ctor assumes ownership,
  // so we don't deallocate the raw pointers
  return mfem::SparseMatrix(row_ptr_raw, col_ind_raw, values_raw, nrows, ncols);
}

// mfem::Vector doesn't have an initializer_list constructor,
// so this is a shorthand to create small vectors in-line
template <typename... types>
mfem::Vector vec(types... args)
{
  mfem::Vector v(sizeof...(types));
  int          count = 0;
  ((v[count++] = args), ...);
  return v;
}

struct Unknown {
  struct Slice {
    Unknown&         x_;
    std::vector<int> ids_;

    Slice(Unknown& x, std::vector<int> ids) : x_(x), ids_(ids) {}

    Slice operator[](int i) { return Slice{x_, {ids_[i]}}; }

    Slice operator[](std::vector<int> i)
    {
      std::vector<int> ids(i.size());
      for (size_t j = 0; j < ids.size(); j++) {
        ids[j] = ids_[i[j]];
      }
      return Slice{x_, ids};
    }

    auto size() { return x_.size_; }
  };

  struct Matvec {
    mfem::SparseMatrix& A_;
    Unknown&            x_;
  };

  Slice operator[](int i) { return Slice{*this, {i}}; }
  Slice operator[](std::vector<int> i) { return Slice{*this, i}; }

  auto size() { return size_; }
  int  size_;
};

// C * x {==, <=, >=} g
struct Constraint {
  enum class Type
  {
    Equality,
    LessThanOrEqual,
    GreaterThanOrEqual
  };

  mfem::SparseMatrix C_;
  mfem::Vector       g_;
  Type               type_;
};

auto operator==(Unknown::Matvec Ax, mfem::Vector& b) { return Constraint{Ax.A_, b, Constraint::Type::Equality}; }
auto operator<=(Unknown::Matvec Ax, mfem::Vector& b) { return Constraint{Ax.A_, b, Constraint::Type::LessThanOrEqual}; }
auto operator>=(Unknown::Matvec Ax, mfem::Vector& b)
{
  return Constraint{Ax.A_, b, Constraint::Type::GreaterThanOrEqual};
}

auto operator==(Unknown::Slice lhs, const mfem::Vector& rhs)
{
  size_t  nrows   = lhs.ids_.size();
  size_t  ncols   = lhs.x_.size();
  size_t  nnz     = nrows;
  int*    row_ptr = new int[nrows + 1];
  int*    col_ind = new int[nnz];
  double* values  = new double[nnz];
  for (size_t i = 0; i < nrows; i++) {
    row_ptr[i] = i;
    col_ind[i] = int(lhs.ids_[i]);
    values[i]  = 1.0;
  }
  row_ptr[nrows] = nnz;
  mfem::SparseMatrix C(row_ptr, col_ind, values, nrows, ncols);

  return Constraint{C, rhs, Constraint::Type::Equality};
}

auto operator==(Unknown::Slice lhs, Unknown::Slice rhs)
{
  size_t  nrows   = lhs.ids_.size();
  size_t  ncols   = lhs.x_.size();
  size_t  nnz     = nrows * 2;
  int*    row_ptr = new int[nrows + 1];
  int*    col_ind = new int[nnz];
  double* values  = new double[nnz];
  for (size_t i = 0; i < nrows; i++) {
    row_ptr[i]         = 2 * i;
    col_ind[2 * i + 0] = int(lhs.ids_[i]);
    values[2 * i + 0]  = 1.0;

    col_ind[2 * i + 1] = int(rhs.ids_[i]);
    values[2 * i + 1]  = -1.0;
  }
  row_ptr[nrows] = nnz;
  mfem::SparseMatrix C(row_ptr, col_ind, values, nrows, ncols);
  mfem::Vector       g(lhs.ids_.size());
  g = 0.0;
  return Constraint{C, g, Constraint::Type::Equality};
}

auto operator*(mfem::SparseMatrix& A, Unknown& x) { return Unknown::Matvec{A, x}; }

struct LinearSystem {
  Constraint              primary_constraint_;
  std::vector<Constraint> auxiliary_constraints_;

  LinearSystem(Constraint p) : primary_constraint_(p), auxiliary_constraints_{} {}

  void append_constraint(Constraint c) { auxiliary_constraints_.push_back(c); }
  void write_out()
  {
    std::ofstream outfile;

    outfile.open("K.mtx");
    primary_constraint_.C_.PrintMM(outfile);
    outfile.close();

    outfile.open("f.vec");
    primary_constraint_.g_.Print(outfile);
    outfile.close();

    for (size_t i = 0; i < auxiliary_constraints_.size(); i++) {
      outfile.open("C" + std::to_string(i) + ".mtx");
      auxiliary_constraints_[i].C_.PrintMM(outfile);
      outfile.close();

      outfile.open("g" + std::to_string(i) + ".vec");
      auxiliary_constraints_[i].g_.Print(outfile);
      outfile.close();
    }
  }
};

int main()
{

  std::string prefix = std::string(SERAC_REPO_DIR) + "/data/matrices/";

  auto compare = [prefix](std::string suffix) {
    std::ofstream outfile(prefix + "output_" + suffix);
    auto K = read_matrix_market(prefix + suffix);
    K.PrintMM(outfile);
    outfile.close();
  };

  compare("bcsstk06.mtx");
  compare("bp_1200.mtx");
  compare("illc1033.mtx");
  compare("plsk1919.mtx");
  compare("young1c.mtx");

#if 0
  auto K = read_sparse_matrix(std::string(SERAC_REPO_DIR) + "/data/K.mtx");
  int  n = K.NumRows();
  auto f = random_vector(n);

  int  m = 10;
  auto C = random_sparse_matrix(n, m, 0.5);
  auto g = random_vector(m);

  Unknown u{n};

  auto sys = LinearSystem(K * u == f);

  // e.g. circuit coupling (2 dof)
  sys.append_constraint(u[{4, 7, 8, 13}] == u[{21, 22, 23, 24}]);

  // e.g. Dirichlet boundary conditions (single dof)
  sys.append_constraint(u[{15, 32, 79, 91}] == vec(1.0, 2.0, 3.0, 4.0));

  // e.g. no-slip condition ({2 or 3}-dof)
  // sys.append_constraint(???);

  // e.g. contact, imcompressible flow/deformation (many dof)
  sys.append_constraint(C * u == g);

  sys.write_out();
#endif
}