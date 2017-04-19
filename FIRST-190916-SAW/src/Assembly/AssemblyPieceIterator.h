#ifndef _ASSEMBLY_PIECE_ITERATOR_
#define _ASSEMBLY_PIECE_ITERATOR_

#include "AssemblyIterator.h"
#include "Piece.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
class AssemblyPieceIterator : public AssemblyIterator<I,T,R> {

private:
  map<I,R> indexToPieceMap;
  typename map<I,R>::iterator indexToPieceMapIter;

public:
  AssemblyPieceIterator( map<I,R> iTPM_ ) :
    indexToPieceMap( iTPM_ ) 
  {
    indexToPieceMapIter = indexToPieceMap.begin();
  };
  ~AssemblyPieceIterator()
  {};

public:
  bool hasNext();
  R next();
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
bool AssemblyPieceIterator<I,T,R>::hasNext(){

  if( indexToPieceMapIter != indexToPieceMap.end() )
    return true;
  
  return false;
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T, class R>
R AssemblyPieceIterator<I,T,R>::next(){

  R piece = indexToPieceMapIter->second;
  ++indexToPieceMapIter;

  return piece;
};

#endif
