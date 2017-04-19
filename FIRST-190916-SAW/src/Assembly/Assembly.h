// Brandon Hespenheide (c) 2007
////////////////////////////////////////////////////////////////////////////////

#ifndef _ASSEMBLY_
#define _ASSEMBLY_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>

#include "Piece.h"
#include "AssemblyPieceIterator.h"
#include "AssemblyStepIterator.h"

using namespace std;

template <class I, class T> class Assembly;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class AssemblyVisitor {

public:
  AssemblyVisitor()
  {};
  virtual ~AssemblyVisitor()
  {};

public:
  virtual void visit( Assembly<I,T> *assembly ){
    cout << "AssemblyVisitor base class not implemented" << endl;
  };
  
};

////////////////////////////////////////////////////////////////////////////////
// Description
////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class Assembly {

protected:
  map< T, Piece<I,T>* > elementToPieceMap;

  // Maps Piece.ID to a piece ptr.
  //////////////////////////////////////////////////////////////////////
  map< PieceID, Piece<I,T>* > id2p;
  typename map< PieceID, Piece<I,T>* >::iterator id2pIter;

  map< I, Piece<I,T>* > indexToPieceMap;

  unsigned int pieceCounter;

public:
  Assembly() :
    pieceCounter(0)
  {};
  ~Assembly()
  {};
    
 public:
  void insertPiece( I index, T *elements );
  unsigned int getTotalPieces() const { return pieceCounter; };

	I getMaxIndex();
	I getMinIndex();
	
  AssemblyPieceIterator<I,T,Piece<I,T>*>* getPieceIterator(){
    return new AssemblyPieceIterator<I,T,Piece<I,T>*>( indexToPieceMap );
  };

  AssemblyStepIterator<I,T,Piece<I,T>*>* getStepIterator(){
    return new AssemblyStepIterator<I,T,Piece<I,T>*>( indexToPieceMap );
  };

  void acceptVisitor( AssemblyVisitor<I,T> *assemblyVisitor );
  void printAssembly();
  void printCutoffs();
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Assembly<I,T>::insertPiece( I index, T *newElements ){

  //cout << "NEW PIECE [" << newElements[0] << " " << newElements[1] << "] : " << index << endl;
  
  // Find all the pieces that overlap with this piece
  //////////////////////////////////////////////////////////////////////  
  map< I, Piece<I,T>* > upperPieces;
  map< T, Piece<I,T>* > lowerPieces;
 
  typename map< T, Piece<I,T>* >::iterator elementToPieceMapIter;

  // Find the piece the new elements belong to, if any. 
  //////////////////////////////////////////////////////////////////////
  for( size_t a = 0; a < 2; a++ ){

    elementToPieceMapIter = elementToPieceMap.find( newElements[a] );
    if( elementToPieceMapIter != elementToPieceMap.end() ) {
      
      Piece<I,T> *currentPiece = elementToPieceMapIter->second;
      
      if( currentPiece->getIndex() > index ){
	//cout << " " << newElements[a] << " is in upper piece " << currentPiece->getID() << endl;
	upperPieces.insert( pair< I, Piece<I,T>* > ( currentPiece->getIndex(), currentPiece ) ); 
      }
      else if( currentPiece->getIndex() < index ){
	
	//cout << " " << newElements[a] << " is in lower piece " << currentPiece->getID() << endl;

	Piece<I,T> *testPiece = currentPiece;
	id2pIter = id2p.find( currentPiece->getLastIndexedUpperPiece() );
	if( id2pIter != id2p.end() &&
	    (id2pIter->second)->getIndex() < index &&
	    (id2pIter->second)->getIndex() > currentPiece->getIndex() ){
	  testPiece = id2pIter->second;
	  //cout << "reset currentPiece " << currentPiece->getID() << " to last indexed piece " << testPiece->getID() << endl;
	}	  
	
	Piece<I,T> *superPiece = testPiece->getSuperPiece();
	
	while( superPiece != NULL &&
	       superPiece->getIndex() < index ) {	       

	  testPiece = superPiece;
	  superPiece = testPiece->getSuperPiece();
	}

	//cout << " Element " << newElements[a] << " is in lower piece " << testPiece->getID() << " " << testPiece->getIndex() << endl; 
	lowerPieces.insert( pair<T, Piece<I,T>* > ( newElements[a], testPiece ) );
	if( testPiece != currentPiece )
	  currentPiece->setLastIndexedUpperPiece( testPiece->getID() );
      }
    }
  }

  // Check for redundancy
  //////////////////////////////////////////////////////////////////////
  if( lowerPieces.size() == 2 &&
      lowerPieces[ newElements[0] ] == lowerPieces[ newElements[1] ] ){
    //cout << "This piece [" << newElements[0] << " " << newElements[1] << "]; " << index << " is redundant at this index value." << endl;
    return;
  }
  
  // Create a new Piece
  //////////////////////////////////////////////////////////////////////
  Piece<I,T> *newPiece = new Piece<I,T>( index );

  if( indexToPieceMap.insert( pair< I, Piece<I,T>* > ( index, newPiece ) ).second ){

    newPiece->addElement( newElements[0] );
    newPiece->addElement( newElements[1] );
    newPiece->setID( ++pieceCounter );
    id2p.insert( pair< int, Piece<I,T>* > ( newPiece->getID(), newPiece ));
  }
  else{

    //cout << "duplicate Index " << index << " [" << newElements[0] << " " << newElements[1] << "]" << endl;
    delete newPiece;
    newPiece = indexToPieceMap[index];

    if( newPiece->hasElement( newElements[0] ) &&
	newPiece->hasElement( newElements[1] ) ){
      //cout << "These elements " <<  newElements[0] << " " << newElements[1] << " already exist at this index value" << endl;
      return;
    }
    
    //cout << "index already exists. Set newPiece to existing piece " << newPiece->getID() << endl;
    newPiece->addElement( newElements[0] );
    //cout << "inserted element " << newElements[0] << " into piece " << newPiece->getID() << endl;
    newPiece->addElement( newElements[1] );
    //cout << "inserted element " << newElements[1] << " into piece " << newPiece->getID() << endl;

    Piece<I,T> *superPiece = newPiece->getSuperPiece();
    if( superPiece != NULL ){
      //cout << " SuperPiece is not NULL" << endl;
      upperPieces.insert( pair< I, Piece<I,T>* > ( superPiece->getIndex(), superPiece ) );
      //cout << " Inserted superPiece " << superPiece->getID() << " into upperPieces map" << endl;
      superPiece->removeSubPiece( newPiece );
    }
  }
  
  // Manage lower pieces
  //////////////////////////////////////////////////////////////////////
  typename map< T, Piece<I,T>* >::iterator lowerPiecesIter = lowerPieces.begin();
  while( lowerPiecesIter != lowerPieces.end() ){

    Piece<I,T> *lowerPiece = lowerPiecesIter->second;

    newPiece->addSubPiece( lowerPiece );
    //cout << " Added piece " << lowerPiece->getID() << " as subpiece of " << newPiece->getID() << endl;

    newPiece->removeElement( lowerPiecesIter->first );
    //cout << " Removed element [" << lowerPiecesIter->first << "] from piece " << newPiece->getID() << endl;
    Piece<I,T> *superPiece = lowerPiece->getSuperPiece();
    if( superPiece != NULL &&
	superPiece != newPiece ){
      //cout << " SuperPiece is not NULL" << endl;
      upperPieces.insert( pair< I, Piece<I,T>* > ( superPiece->getIndex(), superPiece ) );
      //cout << " In remove lower. Inserted superPiece " << superPiece->getID() << " into upperPieces map" << endl;
      superPiece->removeSubPiece( lowerPiece );
      //cout << " Removed " << lowerPiece->getID() << " as subPiece of " << superPiece->getID() << endl;
    }
    //cout << " Adding piece " << newPiece->getID() << " as superPiece of " << lowerPiece->getID() << endl;
    lowerPiece->setSuperPiece( newPiece );
    //cout << " Added piece " << newPiece->getID() << " as superPiece of " << lowerPiece->getID() << endl;

    lowerPiecesIter++;
  }

  //cout << endl;

  // Manage upper pieces
  //////////////////////////////////////////////////////////////////////
  typename map< I, Piece<I,T>* >::iterator upperPieceIter;
  Piece<I,T> *upperPiece = NULL;
  Piece<I,T> *currentPiece = newPiece;

  upperPieces.erase( index );
  
  while( !upperPieces.empty() ){

    upperPieceIter = upperPieces.begin();
    //cout << " Obtained new upperPieceIter " << upperPieces.size() << endl;
    upperPiece = upperPieceIter->second;
    //cout << " Obtained new upperPiece " << upperPiece->getID() << " {" << upperPiece << "}" << endl;

    //cout << " Checking upper Piece " << upperPiece->getID() << " of " << upperPieces.size() << " total upper pieces" << endl;

    // deletion criteria
    //////////////////////////////////////////////////////////////////////
    if( upperPiece->totalSubPieces() == 0 && 
	( upperPiece->totalElements() == 0 ||
	  ( upperPiece->totalElements() == 1 && 
	    ( upperPiece->hasElement( newElements[0] ) || upperPiece->hasElement(newElements[1]) ) ) ) ){
      
      //cout << " Deleting piece " << upperPiece->getID() << endl;

      Piece<I,T> *upperUpperPiece = upperPiece->getSuperPiece();

      if( upperUpperPiece != NULL ){
	//cout << " Checking upperUpperPiece " << upperUpperPiece->getID() << " of piece " << upperPiece->getID() << endl;
	upperUpperPiece->removeSubPiece( upperPiece );
	//cout << " Removed piece " << upperPiece->getID() << " from piece " << upperUpperPiece->getID() << endl;
	upperPieces.insert( pair< I, Piece<I,T>* > ( upperUpperPiece->getIndex(), upperUpperPiece ) );
      }

      indexToPieceMap.erase( upperPieceIter->first );
      //cout << " Removed upperPiece " << upperPiece->getID() << " from indexToPieceMap " << endl;
 
      upperPieces.erase( upperPieceIter );
      //cout << " Removed upperPiece " << upperPiece->getID() << " from upperPieces map" << endl;

      //cout << " Set currentPiece " << currentPiece->getID() << " superPiece to NULL" << endl;
      currentPiece->setSuperPiece( NULL );
      //cout << " Set currentPiece " << currentPiece->getID() << " superPiece to NULL" << endl;
      //cout << " Deleted upperPiece ptr" << endl;

      id2p.erase( upperPiece->getID() );
      //id2p.insert( pair< int, Piece<I,T>* > ( upperPiece->getID(), currentPiece ) );

      //cout << "deleting piece " << upperPiece->getID() << " {" << upperPiece << "}" << endl;
      delete upperPiece;
    }
    else{ // merge
      //cout << " Merging current piece " << currentPiece->getID() << " into upperPiece path " << upperPiece->getID() << endl;

      if( upperPiece->totalElements() != 0 ){
	upperPiece->removeElement( newElements[0] );
	//cout << " Removed element [" << newElements[0] << "] from upperPiece " << upperPiece->getID() << endl;
	upperPiece->removeElement( newElements[1] );
	//cout << " Removed element [" << newElements[1] << "] from upperPiece " << upperPiece->getID() << endl;
      }

      //cout << " Set piece " << upperPiece->getID() << " as superPiece of piece " << currentPiece->getID() << endl;
      currentPiece->setSuperPiece( upperPiece );
      //cout << " Set piece " << upperPiece->getID() << " as superPiece of piece " << currentPiece->getID() << endl;
      upperPiece->addSubPiece( currentPiece );
      //cout << " Set piece " << currentPiece->getID() << " as subPiece of piece " << upperPiece->getID() << endl;
      
      // BMH Stop when there are no other upper pieces. 
      //if( upperPieceIter != upperPieces.end() ){ // if this was the last upper piece, we can stop. 
      //Piece<I,T> *upperUpperPiece = upperPiece->getSuperPiece();
      //if( upperUpperPiece != NULL ){
      //  upperPieces.insert( pair< I, Piece<I,T>* > ( upperUpperPiece->getIndex(), upperUpperPiece ) );
      //  upperUpperPiece->removeSubPiece( upperPiece );
      //cout << " Removed piece " << upperPiece->getID() << " from piece " << upperUpperPiece->getID() << endl;
      //}
      //currentPiece = upperPiece;
      //}

      if( upperPieceIter != upperPieces.end() ){ // if this was the last upper piece, we can stop. 
	
	typename map< I, Piece<I,T>* >::iterator nextUpperPiece = upperPieceIter;
	nextUpperPiece++;

	Piece<I,T> *testPiece = upperPiece;
	id2pIter = id2p.find( upperPiece->getLastIndexedUpperPiece() );
	if( id2pIter != id2p.end() &&
	  (id2pIter->second)->getIndex() < nextUpperPiece->first &&
	  (id2pIter->second)->getIndex() > upperPiece->getIndex() ){
	  testPiece = id2pIter->second;
	  //cout << "reset currentPiece " << currentPiece->getID() << " to last indexed piece " << testPiece->getID() << endl;
	}

	Piece<I,T> *upperUpperPiece = testPiece->getSuperPiece();
	
	while( upperUpperPiece != NULL &&
	       upperUpperPiece->getIndex() < nextUpperPiece->first ){
	  //cout << "checking superPiece " << upperUpperPiece->getID() << " (" << upperUpperPiece->getIndex() << ") of piece " << upperPiece->getID() << endl;
	  testPiece = upperUpperPiece;
	  upperUpperPiece = testPiece->getSuperPiece();
	  //cout << "searching up. upperPiece is " << upperPiece->getID() << " upperUpperPiece is " << upperUpperPiece->getID() << endl;
	}
	
	if( upperUpperPiece != NULL ){
	  
	  upperPieces.insert( pair< I, Piece<I,T>* > ( upperUpperPiece->getIndex(), upperUpperPiece ) );
	  upperUpperPiece->removeSubPiece( testPiece );
	  //cout << " Removed piece " << upperPiece->getID() << " from piece " << upperUpperPiece->getID() << endl;
	}
	currentPiece = testPiece;
      }

      //currentPiece = upperPiece;
      upperPieces.erase( upperPieceIter );
    }
  }

  if( newPiece->hasElement( newElements[0] ) ){
    //cout << " Updating elementMap: " << newElements[0] << " is in piece " << newPiece->getID() << endl;
    elementToPieceMap.erase( newElements[0] );
    elementToPieceMap.insert( pair< T, Piece<I,T>* > ( newElements[0], newPiece ) );
  }
  if( newPiece->hasElement( newElements[1] ) ){
    //cout << " Updating elementMap: " << newElements[1] << " is in piece " << newPiece->getID() << endl;
    elementToPieceMap.erase( newElements[1] );
    elementToPieceMap.insert( pair< T, Piece<I,T>* > ( newElements[1], newPiece ) );
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Assembly<I,T>::printAssembly(){
  int counter = 1;
  typename map< I, Piece<I,T>* >::iterator indexToPieceMapIter = indexToPieceMap.begin();
  while( indexToPieceMapIter != indexToPieceMap.end() ){
    cout << "Step " << counter++ << endl;
    indexToPieceMapIter->second->printPiece();
    cout << endl;
    indexToPieceMapIter++;
  }

}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Assembly<I,T>::printCutoffs(){

  int multiPieceIndex = 0;
  Piece<I,T> *testPiece;

  ofstream newCutoffs;
  newCutoffs.open( "assemblyCutoffs.txt" );

  typename map< I, Piece<I,T>* >::iterator indexToPieceMapIter = indexToPieceMap.begin();
  while( indexToPieceMapIter != indexToPieceMap.end() ){
    
    newCutoffs << showpoint << setw(16) << setprecision(16) << indexToPieceMapIter->first << endl;

    testPiece = indexToPieceMapIter->second;
    //if( testPiece->totalSubPieces() + testPiece->totalElements() -2 > 0 ){
    multiPieceIndex += testPiece->totalSubPieces() + testPiece->totalElements() -2;     
    //cout << " found duplicate index piece? " << testPiece->getID() << endl;
    //testPiece->printPiece();
    //}

    indexToPieceMapIter++;
  }

  //cout << "Total indecies found = " << indexToPieceMap.size() << endl;
  //cout << "Duplicate Index Count = " << multiPieceIndex << endl;

};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
I Assembly<I,T>::getMaxIndex() {
	typename map< I, Piece<I,T>* >::reverse_iterator indexToPieceMapIter = indexToPieceMap.rbegin();

	I maxIndex = (*indexToPieceMapIter).first;
	
	return maxIndex;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
I Assembly<I,T>::getMinIndex() {
	typename map< I, Piece<I,T>* >::iterator indexToPieceMapIter = indexToPieceMap.begin();
	
	I minIndex = (*indexToPieceMapIter).first;
	
	return minIndex;
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void Assembly<I,T>::acceptVisitor( AssemblyVisitor<I,T> *assemblyVisitor ) {
  
  assemblyVisitor->visit( this );
};

#endif
