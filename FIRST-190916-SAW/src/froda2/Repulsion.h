#ifndef REPULSION_H_
#define REPULSION_H_

class Repulsion
{
public:
	Repulsion() {}
	virtual ~Repulsion() {}
  virtual double getInteractionCutoff( int p1, int p2 ) const = 0;
  double operator()( int p1, int p2 ) const { return getInteractionCutoff( p1, p2 ); }
  virtual double getMaxInteractionCutoff( int p ) const = 0;
};

class DefaultRepulsion : public Repulsion
{
public:
  DefaultRepulsion() : Repulsion() {}
  ~DefaultRepulsion() {}
  double getInteractionCutoff( int p1, int p2 ) const { return 1.0; }
  double getMaxInteractionCutoff( int p ) const { return 1.0; }
};
#endif /*REPULSION_H_*/
