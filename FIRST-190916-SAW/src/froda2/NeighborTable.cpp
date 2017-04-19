#include "NeighborTable.h"
#include <map>
#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

void NeighborTable::insert( int index1, int index2 ) {
  //inserts and preserves sorted lists
  //insertion is symmetric
  
  vector<int>::iterator insertionPoint;
  vector<int>* neighborsOfP;

  neighborsOfP = &table1[index1];
  insertionPoint = lower_bound( neighborsOfP->begin(), neighborsOfP->end(), index2 );
  if ( insertionPoint == neighborsOfP->end() || *insertionPoint != index2 )
    neighborsOfP->insert( insertionPoint, index2 );

  //repeat
  neighborsOfP = &table1[index2];
  insertionPoint = lower_bound( neighborsOfP->begin(), neighborsOfP->end(), index1 );
  if ( insertionPoint == neighborsOfP->end() || *insertionPoint != index1 )
    neighborsOfP->insert( insertionPoint, index1 );
}

NeighborTable::NeighborTable( int Npoints ) :
  table2_(NULL),
  table3_(NULL),
  table1( Npoints )
{
}

void NeighborTable::insert( const NeighborTable& otherNeighborTable ) {
  if ( size() != otherNeighborTable.size() ) {
    cout << "Cannot combine tables.  Not the same size" << endl;
    exit(0);
  }
  
  int nP = otherNeighborTable.size();
  for ( int p1 = 0; p1 < nP; p1++ ) {
    for ( size_t j = 0; j < otherNeighborTable[p1].size(); j++ ) {
      int p2 = otherNeighborTable[p1][j];
      //we know that the insert function works symmetrically,
      //so we only need to insert the pair once.
      if ( p1 < p2 ) insert( p1, p2 );
    }
  }
  
}   

ostream& operator<< (ostream& os, const vector< set<int> > &vecsetint ) {
  //const set<int> *pointlist;
  for ( size_t p=0; p<vecsetint.size(); p++ ) {
    os << p << " | ";
    set<int>::const_iterator i;
    for ( i = vecsetint[p].begin(); i != vecsetint[p].end(); i++ ) {
      os << *i << " ";
    }
    os << '\n';
  }
  return os;
}

NeighborTable::~NeighborTable()
{
}

void NeighborTable::setupSecondNeighborTable() {
  int N = size();
  for ( int p=0; p<N; p++ ) {
    const vector<int> *allFirstNeighborsOfP = &table1[p];
    set<int> *allSecondNeighborsOfP = &(*table2_)[p];

    vector<int>::const_iterator n1;
    for ( n1 = allFirstNeighborsOfP->begin(); n1 != allFirstNeighborsOfP->end(); n1++ ) {
      // we can visit all of the first neighbors of p,
      // no need to exclude any.
      // We assume that p is not found in its own first neighbor list.
      const vector<int> *secondNeighbors = &table1[*n1];

      vector<int>::const_iterator n2;
      for ( n2 = secondNeighbors->begin(); n2 != secondNeighbors->end(); n2++ ) {
        // of the second neighbors of p, we must exclude any that happen to
        // be p or a first neighbor of p.
        if ( *n2 == p ) continue;
        if ( find( allFirstNeighborsOfP->begin(), allFirstNeighborsOfP->end(), *n2 ) 
              != allFirstNeighborsOfP->end() ) continue;
        //if ( allFirstNeighborsOfP->find( *n2 ) != allFirstNeighborsOfP->end() ) continue;
        
        // because the second neighbor table is a set, inserting a duplicate
        // will have no effect
        allSecondNeighborsOfP->insert( *n2 );
      }
    }
  }      
}

void NeighborTable::setupThirdNeighborTable() {
  int N = size();
  for ( int p=0; p<N; p++ ) {
    const vector<int> *allFirstNeighborsOfP = &table1[p];
    const set<int> *allSecondNeighborsOfP = &(*table2_)[p];
    set<int> *allThirdNeighborsOfP = &(*table3_)[p];
    
    set<int>::const_iterator n2;
    for ( n2 = allSecondNeighborsOfP->begin(); n2 != allSecondNeighborsOfP->end(); n2++ ) {
      // we can visit all of the second neighbors of p,
      // no need to exclude any.  We know that this second neighbor
      // list only has the correct second neighbors.      
      const vector<int> *thirdNeighbors = &table1[*n2];
      vector<int>::const_iterator n3;
      for ( n3 = thirdNeighbors->begin(); n3 != thirdNeighbors->end(); n3++ ) {
        // We must exclude any third neighbors that happen to be p,
        // or a first or second neighbor of p.
        if ( *n3 == p ||
             find( allFirstNeighborsOfP->begin(), allFirstNeighborsOfP->end(), *n3 ) 
               != allFirstNeighborsOfP->end() ||
             //allFirstNeighborsOfP->find( *n3 ) != allFirstNeighborsOfP->end() ||
             allSecondNeighborsOfP->find( *n3 ) != allSecondNeighborsOfP->end() )
          continue;
          
        allThirdNeighborsOfP->insert( *n3 );
      }
    }  
  }  
}

void NeighborTable::commit() {

  table2_ = new vector< set<int> >( size() );
  table3_ = new vector< set<int> >( size() );
  setupSecondNeighborTable();
  setupThirdNeighborTable();

  int nP = size();
  table3_onlyPairsILessThanJ.resize(nP);
  for ( int i=0; i<nP; i++ ) {
    set<int>::const_iterator firstJGreaterThanI = (*table3_)[i].upper_bound(i);
    set<int>::const_iterator end = (*table3_)[i].end();
    set<int>::const_iterator iter;
    int count;
    for ( count = 0, iter = firstJGreaterThanI;
          iter != end;
          iter++, count++ ) {}
    table3_onlyPairsILessThanJ[i].resize(count);
    copy( firstJGreaterThanI, end, table3_onlyPairsILessThanJ[i].begin() );
  }

  delete table2_;
  delete table3_;

  /*
  int s=0;
  size_t smax = 0;
  for ( int p=0; p<nP; p++ ) {
    s += table3_onlyPairsILessThanJ[p].size();
    smax = ( table3_onlyPairsILessThanJ[p].size() > smax ) ? table3_onlyPairsILessThanJ[p].size() : smax;
  }
  double s_mean = static_cast<double>(s) / nP;
  cout << "Mean number of 3rd neighbors j>i " << s_mean << endl;
  cout << "Max number of 3rd neighbors j>i " << smax << endl;
  */
  
}
