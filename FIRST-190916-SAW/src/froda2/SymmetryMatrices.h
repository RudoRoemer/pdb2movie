#ifndef SYMMETRY_MATRICES_H_
#define SYMMETRY_MATRICES_H_

#include <vector>
#include "SymReplica.h"

class SymmetryMatrices {
public:
  SymmetryMatrices(vector <Matrix> *matrices_);
  ~SymmetryMatrices();
  const Matrix & getMatrix(int n) const;
  const Matrix & getInverse(int n) const;
  const Matrix & getProduct(int n, int m) const;
  const Vec3 & getOffsetPair(int n, int m) const;
  size_t size() const;
private:
  vector <Matrix> *matrices;
  vector <Matrix> *inverses;
  vector < vector <Matrix> > *products;
  // contains all pairwise products of matrices and inverses
  vector < vector <Vec3> > *offsetPairs;
  // contains all pairwise differences of offset vectors of matrices
  vector <Vec3> averagePos;
  // contains symmetry-averaged positions for all atoms
  Matrix mult(const Matrix &m1, const Matrix &m2);
  Vec3 subtractOffsets(const Matrix &m1, const Matrix &m2);
  double determinant(const Matrix &m);
  void showMatrix(const Matrix &m);
  Matrix invert(const Matrix &m);
};

#endif /* SYMMETRY_MATRICES_H_ */
