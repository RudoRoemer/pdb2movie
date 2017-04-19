#include "RMSD.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include "Vec3.h"

using namespace std;

RMSD::RMSD( const vector<Vec3> *points1_, const vector<Vec3> *points2_ ) :
  points1(points1_),
  points2(points2_)
{
}

RMSD::~RMSD()
{
}

double RMSD::calc() const {
  int nP = points1->size();
  if ( nP != static_cast<int>(points2->size()) ) {
    cout << "Error: RMSD point sets have different numbers of points" << endl;
    exit(0);
  }
  double dist2sum = 0.0;
  for ( int p=0; p<nP; p++ ) {
    dist2sum += (*points1)[p].dist2( (*points2)[p] );
  }
  return sqrt(dist2sum/nP);
}
