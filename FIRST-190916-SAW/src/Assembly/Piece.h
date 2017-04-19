#ifndef _PIECE_
#define _PIECE_

// Brandon Hespenheide (c) 2007
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>

using namespace std;

typedef unsigned int PieceID;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class Piece {
  
protected:
  set<T> elements;

	// cache for allElements 
	// WARNING - there is currently no checking to verify that the sub-elements
	//           have not changed
	// FIXME   - implement some form of checking to correct that
	mutable set<T> allElements;
	
  Piece<I,T> *superPiece;
  set< Piece<I,T>* > subPieces;
  I index;
  PieceID id; 
  PieceID lastIndexedUpperPiece;

public:
  Piece( I i, set<T> initialElements ) :
    elements( initialElements ), 
    superPiece(NULL),
    index(i),
    id(0),
    lastIndexedUpperPiece(0)
  {};
  Piece( I i ) :
    superPiece(NULL),
    index(i),
    id(0),
    lastIndexedUpperPiece(0)
  {};
  ~Piece(){};

public:

  void setIndex( I i );
  I getIndex();

  void setID( PieceID i );
  PieceID getID();

  void setLastIndexedUpperPiece( int i );
  int getLastIndexedUpperPiece();

  void setSuperPiece( Piece<I,T> *superP );
  Piece<I,T>* getSuperPiece();

  void addSubPiece( Piece<I,T> *subP );
  void removeSubPiece( Piece<I,T> *subP );
  size_t totalSubPieces();

  void addElement( T e );
  set<T> getElements() const;
  set<T> getElementsAndSubElements() const;
  void removeElement( T e );
  bool hasElement( T e );
  size_t totalElements();

  void printPiece();
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::setIndex( I i ){

  index = i;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
I Piece<I,T>::getIndex(){ 

  return index; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::setID( PieceID i ){ 
  
  id = i; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
PieceID Piece<I,T>::getID(){ 

  return id; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::setLastIndexedUpperPiece( int i ){ 
  
  lastIndexedUpperPiece = i;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
int Piece<I,T>::getLastIndexedUpperPiece(){ 
  
  return lastIndexedUpperPiece; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
Piece<I,T>* Piece<I,T>::getSuperPiece(){ 

  return superPiece; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::setSuperPiece( Piece<I,T> *superP ){ 

  if( superP != NULL &&
      superP->getID() == id ){
    cout << "attempting to set superPiece to self on piece " << id << endl;
    exit(1);
  }
  superPiece = superP;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
size_t Piece<I,T>::totalSubPieces(){ 
  
  return subPieces.size(); 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::addSubPiece( Piece<I,T> *subP ){ 

  if( subP->getID() == id ){
    cout << "attempting to set subPiece to self on piece " << id << endl;
    exit(1);
  }
  subPieces.insert( subP ); 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::removeSubPiece( Piece<I,T> *subP ){

  subPieces.erase( subP );

  /*
  cout << "checking " << subPieces.size() << " subPieces" << endl;
  typename set< Piece<I,T>* >::iterator subPiecesIter = subPieces.begin();
  while( subPiecesIter != subPieces.end() ){
  
    cout << subP << "  " << *subPiecesIter << " id " << (*subPiecesIter)->getID() << endl;

    if( *subPiecesIter == subP ){
      subPieces.erase( subPiecesIter );
    }
    else{
    subPiecesIter++;
    }
  }
  */
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
size_t Piece<I,T>::totalElements(){ 

  return elements.size();
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
set<T> Piece<I,T>::getElements() const { 

  return elements; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
set<T> Piece<I,T>::getElementsAndSubElements() const {
  
	if (allElements.size() == 0) { // FIXME - add a check to recompute allElements 
																 // when necessary (ie - if a subPiece has 
																 // been added or removed)
		// initialize with elements in this piece
		allElements.insert(elements.begin(), elements.end());
		
		typename set<Piece<I,T>*>::iterator subPiecesIter = subPieces.begin();
		while( subPiecesIter != subPieces.end() ){
			set<T> subPieceElements = (*subPiecesIter)->getElementsAndSubElements();
			allElements.insert( subPieceElements.begin(), subPieceElements.end() );      
			subPiecesIter++;
		}
	}
	
	return allElements; 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::addElement( T e ){ 

  elements.insert(e); 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::removeElement( T e ){ 

  elements.erase( e ); 
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
bool Piece<I,T>::hasElement( T e ){

  if( elements.find(e) != elements.end() )
    return true;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Piece<I,T>::printPiece(){
  
  typename set<T>::iterator print = elements.begin();
  cout << index << " [P ";
  while( print != elements.end() ){
    cout << *print << " ";
    print++;
  }
  
  typename set< Piece<I,T>* >::iterator subPiecesIter = subPieces.begin();
  while( subPiecesIter != subPieces.end() ){
    cout << " [S ";
    (*subPiecesIter)->printPiece();
    subPiecesIter++;
    cout << "]";
  }
  cout << "]";
}

#endif
