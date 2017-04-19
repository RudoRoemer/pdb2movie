#ifndef ROPESASSIGNMENT_H_
#define ROPESASSIGNMENT_H_

class Rope {
public:
  Rope() {}
  Rope( int p1_, int p2_, double length_ ) :
    p1(p1_),
    p2(p2_),
    length(length_) {}
  Rope& operator=( const Rope& other ) {
    p1 = other.p1;
    p2 = other.p2;
    length = other.length;
    return *this;
  }
  int p1;
  int p2;
  double length;
};

#include <vector>

class RopesAssignment
{
public:
	RopesAssignment() {}
	virtual ~RopesAssignment() {}
  
  void insert( const Rope& rope ) {
    ropes.push_back( rope );
  }

  void insert( int p1, int p2, double length ) {
    ropes.push_back( Rope( p1, p2, length ) );
  }
  
  void insert( const RopesAssignment& otherRopesAssignment ) {
    ropes.insert( ropes.end(),
                  otherRopesAssignment.begin(),
                  otherRopesAssignment.end() );
  }
  
  typedef std::vector<Rope>::iterator iterator;
  typedef std::vector<Rope>::const_iterator const_iterator;

  iterator begin() { return ropes.begin(); }
  const_iterator begin() const { return ropes.begin(); }
  iterator end() { return ropes.end(); }
  const_iterator end() const { return ropes.end(); }
  size_t size() const { return ropes.size(); }
  Rope& operator[]( size_t i ) { return ropes[i]; }
  const Rope& operator[]( size_t i ) const { return ropes[i]; }
  
private:
  std::vector<Rope> ropes;
};

#endif /*ROPESASSIGNMENT_H_*/
