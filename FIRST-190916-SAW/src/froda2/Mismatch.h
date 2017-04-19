#ifndef MISMATCH_H_
#define MISMATCH_H_

#include <list>
#include "Observable.h"

class MismatchTerm
{
public:
  MismatchTerm() : verbose( false ) {}
  virtual ~MismatchTerm() {}  
  virtual double mismatch() = 0;
  
  void setVerbose( bool val ) {
    verbose = val;
  }
protected:
  bool verbose;  
};

class Mismatch : public Observer {
public:
  Mismatch();
  virtual ~Mismatch();
  void addTerm( MismatchTerm* term ) { mismatchTerms.push_back(term); }
  double calc() {
    if ( !isUpdated ) update();
    return maxMismatch;
  }
  double operator()() { return calc(); }
  void setVerbose( bool val ) {
    for ( std::list<MismatchTerm*>::iterator it = mismatchTerms.begin();
          it != mismatchTerms.end();
          it++ ) {
      (*it)->setVerbose( val );  
    }
    isUpdated = false;
  }
  void receiveNotification( Observable *observable );
  void update();
private:
  std::list<MismatchTerm*> mismatchTerms;
  double maxMismatch;
  bool isUpdated;
};

#endif /*MISMATCH_H_*/
