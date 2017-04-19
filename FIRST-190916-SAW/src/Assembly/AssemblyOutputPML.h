// Brandon Hespenheide (c) 2007
////////////////////////////////////////////////////////////////////////////////

#ifndef _ASSEMBLY_OUTPUT_PML_
#define _ASSEMBLY_OUTPUT_PML_

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Assembly.h"
#include "AssemblyStepIterator.h"
#include "MolFramework.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class AssemblyOutputPML : public AssemblyVisitor<I,T> {
	
private:
	MolFramework *molFramework;
		
public:
		AssemblyOutputPML(MolFramework *_molFramework ) :
    molFramework( _molFramework ) {	};
  ~AssemblyOutputPML()
  {};
	
public:
		virtual void visit( Assembly<I,T> *assembly );
	
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void AssemblyOutputPML<I,T>::visit( Assembly<I,T> *assembly ){
	
	setupDimensions(assembly);
	
  cout << "AssemblyOutputPML currently being implemented." << endl;
  cout << "molFramework pointer has been set. Total sites: " << molFramework->total_sites << endl;
  
	
};

#endif
