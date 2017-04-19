#ifndef COVALENTNEIGHBORTABLEBUILDER_H_
#define COVALENTNEIGHBORTABLEBUILDER_H_

#include "NeighborTable.h"

class CovalentNeighborTableBuilder
{
public:
	CovalentNeighborTableBuilder() {}
	virtual ~CovalentNeighborTableBuilder() {}
  virtual NeighborTable *getNeighborTable() = 0; 
};

#endif /*COVALENTNEIGHBORTABLEBUILDER_H_*/
