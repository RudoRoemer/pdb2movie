// SymmetryMatrices.cpp
// Stores matrices, inverses, and their products
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include "SymmetryMatrices.h"

using namespace std;

SymmetryMatrices::SymmetryMatrices(vector <Matrix> *matrices_) {
  matrices = matrices_;
  inverses = new vector <Matrix>;
  products = new vector < vector <Matrix> >;
  offsetPairs = new vector < vector <Vec3> >;
  for (size_t i = 0; i < matrices->size(); i++) {
    inverses->push_back(invert(matrices->at(i)));
  }
  products->resize(matrices->size());
  offsetPairs->resize(matrices->size());
  for (size_t i = 0; i < matrices->size(); i++) {
    products->at(i).resize(matrices->size());
    offsetPairs->at(i).resize(matrices->size());
    for (size_t j = 0; j < matrices->size(); j++) {
      products->at(i).at(j) = mult(matrices->at(i),inverses->at(j));
      offsetPairs->at(i).at(j) = subtractOffsets(matrices->at(i),matrices->at(j));
    }
  }
}

SymmetryMatrices::~SymmetryMatrices() {
  delete inverses;
  delete products;
}


Matrix SymmetryMatrices::mult(const Matrix &m1, const Matrix &m2) {
  // For calculating gradients I only need the rotational parts of the
  // matrices, not the constant offset vectors
  Matrix result; // make this 3x3
  result.resize(3);
  for (int i = 0; i < 3; i++) {
    result[i].resize(3);
    for (int j = 0; j < 3; j++) {
      result[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        result[i][j] += m1[i][k]*m2[k][j];
      }
    }
  }
  return result;
}

Vec3 SymmetryMatrices::subtractOffsets(const Matrix &m1, const Matrix &m2) {
  // if o(n) is the vector offset and m(n) is the rotation matrix, I want
  // o(m1) - m(m1)o(m2)
  Vec3 result;
  result.x = m1[0][3];
  result.y = m1[1][3];
  result.z = m1[2][3];
  for (int k = 0; k < 3; k++) {
    result.x -= m1[k][0]*m2[0][3];
    result.y -= m1[k][1]*m2[1][3];
    result.z -= m1[k][2]*m2[2][3];
  }
  return result;
}


double SymmetryMatrices::determinant(const Matrix &m) {
  // Assumes it's being given a 3x3 matrix
  double result = 0;
  result += m[0][0]*m[1][1]*m[2][2];
  result -= m[0][0]*m[1][2]*m[2][1];
  result -= m[0][1]*m[1][0]*m[2][2];
  result += m[0][1]*m[1][2]*m[2][0];
  result += m[0][2]*m[1][0]*m[2][1];
  result -= m[0][2]*m[1][1]*m[2][0];
  // there's probably a more elegant way to do that
  return result;
}

void SymmetryMatrices::showMatrix(const Matrix &m) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << setw(8) << m[i][j];
    }
    cout << endl;
  }
  return;
}

Matrix SymmetryMatrices::invert(const Matrix &m) {
  double det = determinant(m);
  if (det == 0) {
    cerr << "ERROR: Symmetry matrix has no inverse!\n";
    showMatrix(m);
    exit(EXIT_FAILURE);
  }
  Matrix result;
  result.resize(3);
  for (int i = 0; i < 3; i++) {
    result[i].resize(3);
  }
  result[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  result[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
  result[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
  result[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  result[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
  result[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
  result[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  result[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
  result[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      result[i][j] /= det;
    }
  }
  return result;
}

const Matrix & SymmetryMatrices::getMatrix(int n) const {
  return matrices->at(n);
}

const Matrix & SymmetryMatrices::getInverse(int n) const {
  return inverses->at(n);
}

const Matrix & SymmetryMatrices::getProduct(int n, int m) const {
  return products->at(n)[m];
}

const Vec3 & SymmetryMatrices::getOffsetPair(int n, int m) const {
  return offsetPairs->at(n)[m];
}

size_t SymmetryMatrices::size() const {
  return matrices->size();
}

