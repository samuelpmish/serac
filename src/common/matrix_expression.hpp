// Copyright (c) 2019-2020, Lawrence Livermore National Security, LLC and
// other Serac Project Developers. See the top-level LICENSE file for
// details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * @file matrix_expression.hpp
 *
 * @brief A set of template classes used to represent the evaluation of unary and binary operations
 * on matrices
 */

#ifndef MATRIX_EXPRESSION
#define MATRIX_EXPRESSION

#include "mfem.hpp"

namespace serac {
/**
 * @brief A base class representing a matrix expression
 * @tparam T The base vector type, e.g., mfem::DenseMatrix, or another MatrixExpr
 * @note This class should never be used directly
 */
template <typename T>
class MatrixExpr {
public:
  /**
   * @brief Returns the fully evaluated value for the matrix
   * expression at index @p i, @p j
   * @param i The row index to evaluate at
   * @param j The col index to evaluate at
   */
  double operator()(size_t i, size_t j) const { return asDerived()(i, j); }

  /**
   * @brief Returns the width (cols) of the matrix expression
   */
  size_t Width() const { return asDerived().Width(); }

  /**
   * @brief Returns the height (rows) of the matrix expression
   */
  size_t Height() const { return asDerived().Height(); }

  /**
   * @brief Implicit conversion operator for fully evaluating
   * a matrix expression into an actual mfem::DenseMatrix
   * @return The fully evaluated matrix
   */
  operator mfem::DenseMatrix() const
  {
    mfem::DenseMatrix result(Height(), Width());
    for (size_t i = 0; i < Height(); i++) {
      for (size_t j = 0; j < Width(); j++) {
        result(i, j) = (*this)(i, j);
      }
    }
    return result;
  }

  /**
   * @brief Performs a compile-time downcast to the derived object
   * @return The derived object
   * @see Curiously Recurring Template Pattern
   */
  const T& asDerived() const { return static_cast<const T&>(*this); }

  /**
   * @brief Performs a compile-time downcast to the derived object
   * @return The derived object
   * @see Curiously Recurring Template Pattern
   */
  T& asDerived() { return static_cast<T&>(*this); }
};

/**
 * @brief Fully evaluates a matrix expression into an actual mfem::DenseMatrix
 * @return The fully evaluated matrix
 * @see MatrixExpr::operator mfem::DenseMatrix
 */
template <typename T>
mfem::DenseMatrix evaluate(const MatrixExpr<T>& expr)
{
  return expr;
}

}  // namespace serac

#endif
