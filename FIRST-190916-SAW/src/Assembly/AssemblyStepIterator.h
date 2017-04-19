#ifndef _ASSEMBLY_STEP_ITERATOR_
#define _ASSEMBLY_STEP_ITERATOR_

#include <vector>
#include <cassert>

#include "AssemblyIterator.h"
#include "Piece.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
class AssemblyStepIterator : public AssemblyIterator<I,T,vector<R> > {

private:
  map<I,R> indexToPieceMap;
  typename map<I,R>::iterator indexToPieceMapIter;
  vector<R> piecesInCurrentStep;

public:
  AssemblyStepIterator( map<I,R> iTPM_ ) :
    indexToPieceMap( iTPM_ )
  {
    indexToPieceMapIter = indexToPieceMap.begin();
  };
  ~AssemblyStepIterator()
  {};

public:
  bool hasNext();
  vector<R> next();

  // Member fxn's unique to this class
  //////////////////////////////////////////////////////////////////////
  vector<R> gotoStepWithIndex( I index );
  I getCurrentStepIndex();
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
bool AssemblyStepIterator<I,T,R>::hasNext(){

  if( indexToPieceMapIter != indexToPieceMap.end() )
    return true;
  
  return false;
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
vector<R> AssemblyStepIterator<I,T,R>::next(){

  // The next peice in the list is automatically added to the 'pieces'
  // in the current step
  //////////////////////////////////////////////////////////////////////
  Piece<I,T> *nextPiece = indexToPieceMapIter->second;
  piecesInCurrentStep.push_back( nextPiece );
  ++indexToPieceMapIter;

  // We need to remove any subpieces of the piece we just added
  //////////////////////////////////////////////////////////////////////
  typename vector<R>::iterator piecesInCurrentStepIter = piecesInCurrentStep.begin();
  while( piecesInCurrentStepIter != piecesInCurrentStep.end() ){

    if( (*piecesInCurrentStepIter)->getSuperPiece() == nextPiece ){
      piecesInCurrentStep.erase( piecesInCurrentStepIter );
    }
    else {
      ++piecesInCurrentStepIter;
    }
  }

  return piecesInCurrentStep;
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
vector<R> AssemblyStepIterator<I,T,R>::gotoStepWithIndex( I index ){

  // Since we don't know where this step is relative to the iterator, we
  // need to clear the pieces in the 'last' step and start over. 
  //////////////////////////////////////////////////////////////////////
  piecesInCurrentStep.clear();

  // BMH implement this. Must return something and use internal iterator in 
  // order to keep next() and hasNext() fxns working properly. 

}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
I AssemblyStepIterator<I,T,R>::getCurrentStepIndex(){

  // You have to call the fxn 'next()' to create a currentStep before you
  // can get the currentStepIndex...
  //////////////////////////////////////////////////////////////////////
  assert( !piecesInCurrentStep.empty() );
  
  R lastPieceAddedToStep = piecesInCurrentStep.at( piecesInCurrentStep.size() - 1);
  return lastPieceAddedToStep->getIndex();
}

#endif
